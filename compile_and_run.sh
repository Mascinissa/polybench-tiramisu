#!/bin/bash
# Standalone build & measure for one PolyBench-Tiramisu benchmark.
#
#   ./compile_and_run.sh <benchmark> <SIZE> [n_measurements]
#
#   <benchmark>       e.g. gemm, 2mm, floyd_warshall, jacobi1d, heat3d
#   <SIZE>            MINI | SMALL | MEDIUM | LARGE | XLARGE
#   [n_measurements]  default 5 (PolyBench's time_benchmark.sh count)
#
# Measurement methodology is PolyBench/C 4.2.1's: page-aligned
# allocation, deterministic init_array, cache flush before the timer, one
# kernel invocation per fresh process. Times are printed in milliseconds,
# followed by the PolyBench-normalized time (drop min & max, mean of the
# rest) when n_measurements >= 3.
#
# Environment:
#   TIRAMISU_ROOT       Tiramisu installation (required unless TPB_* set).
#                       Both the classic source layout (3rdParty/...) and
#                       the install layout (install/include, install/lib*)
#                       are detected.
#   TPB_INCLUDES        colon-separated include dirs (overrides detection)
#   TPB_LIBS            colon-separated lib dirs (overrides detection)
#   TPB_EXTRA_INCLUDES  colon-separated dirs appended to the include path
#   TPB_EXTRA_LIBS      colon-separated dirs appended to the lib path
#   CXX                 compiler (default g++)
#   TPB_DUMP=1          dump live-out arrays (PolyBench format) to
#                       <function>.dump for correctness checking
#   TPB_KEEP=1          keep build artifacts

set -e

if [ $# -lt 2 ]; then
    echo "Usage: ./compile_and_run.sh <benchmark> <SIZE> [n_measurements]"
    echo "Example: ./compile_and_run.sh gemm MEDIUM 5"
    exit 1
fi

KERNEL=$1
SIZE=$2
RUNS=${3:-5}
[ "${SIZE}" = "EXTRALARGE" ] && SIZE=XLARGE

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
FUNC="function_${KERNEL}_${SIZE}"

GENERATOR=$(find "${REPO_ROOT}" -path "${REPO_ROOT}/.git" -prune -o -name "${FUNC}_generator.cpp" -print | head -n 1)
if [ -z "${GENERATOR}" ]; then
    echo "Benchmark variant not found: ${FUNC}"
    echo "Available benchmarks:"
    find "${REPO_ROOT}" -path "${REPO_ROOT}/.git" -prune -o -name "function_*_MINI_generator.cpp" -print \
        | sed 's/.*function_\(.*\)_MINI_generator.cpp/  \1/' | sort
    exit 1
fi
KERNEL_FOLDER="$(dirname "${GENERATOR}")"

# ---------------------------------------------------------------------------
# Locate Tiramisu (includes + libs)
# ---------------------------------------------------------------------------
if [ -n "${TPB_INCLUDES}" ] && [ -n "${TPB_LIBS}" ]; then
    INCLUDE_DIRS="${TPB_INCLUDES}"
    LIB_DIRS="${TPB_LIBS}"
elif [ -n "${TIRAMISU_ROOT}" ] && [ -d "${TIRAMISU_ROOT}/install/include" ]; then
    # install layout
    INCLUDE_DIRS="${TIRAMISU_ROOT}/install/include"
    LIB_DIRS="${TIRAMISU_ROOT}/install/lib:${TIRAMISU_ROOT}/install/lib64"
elif [ -n "${TIRAMISU_ROOT}" ] && [ -d "${TIRAMISU_ROOT}/3rdParty" ]; then
    # classic source layout; tolerate both Halide dir arrangements
    INCLUDE_DIRS="${TIRAMISU_ROOT}/include"
    LIB_DIRS="${TIRAMISU_ROOT}/build"
    for d in "${TIRAMISU_ROOT}/3rdParty/Halide/include" \
             "${TIRAMISU_ROOT}/3rdParty/Halide/install/include" \
             "${TIRAMISU_ROOT}/3rdParty/isl/include" \
             "${TIRAMISU_ROOT}/3rdParty/isl/build/include"; do
        [ -d "$d" ] && INCLUDE_DIRS="${INCLUDE_DIRS}:$d"
    done
    for d in "${TIRAMISU_ROOT}/3rdParty/Halide/lib" \
             "${TIRAMISU_ROOT}/3rdParty/Halide/lib64" \
             "${TIRAMISU_ROOT}/3rdParty/Halide/install/lib" \
             "${TIRAMISU_ROOT}/3rdParty/Halide/install/lib64" \
             "${TIRAMISU_ROOT}/3rdParty/isl/build/lib"; do
        [ -d "$d" ] && LIB_DIRS="${LIB_DIRS}:$d"
    done
else
    echo "Tiramisu not found. Either:"
    echo "  export TIRAMISU_ROOT=<path/to/tiramisu>"
    echo "or set TPB_INCLUDES / TPB_LIBS explicitly."
    exit 1
fi
[ -n "${TPB_EXTRA_INCLUDES}" ] && INCLUDE_DIRS="${INCLUDE_DIRS}:${TPB_EXTRA_INCLUDES}"
[ -n "${TPB_EXTRA_LIBS}" ] && LIB_DIRS="${LIB_DIRS}:${TPB_EXTRA_LIBS}"

CXX=${CXX:-g++}
CXXFLAGS="-std=c++17 -O3 -fno-rtti"

INCLUDES="-I${REPO_ROOT}/utilities -I${KERNEL_FOLDER}"
for d in $(echo "${INCLUDE_DIRS}" | tr ':' ' '); do INCLUDES="${INCLUDES} -I$d"; done
LIB_FLAGS=""
for d in $(echo "${LIB_DIRS}" | tr ':' ' '); do LIB_FLAGS="${LIB_FLAGS} -L$d"; done
export LD_LIBRARY_PATH="${LIB_DIRS}:${LD_LIBRARY_PATH}"

cd "${KERNEL_FOLDER}"

# ---------------------------------------------------------------------------
# 1. Compile & run the Tiramisu generator -> ${FUNC}.o
# ---------------------------------------------------------------------------
echo "[1/3] Compiling and running the Tiramisu generator (${FUNC})"
${CXX} ${CXXFLAGS} ${INCLUDES} "${FUNC}_generator.cpp" ${LIB_FLAGS} \
    -ltiramisu -lHalide -lisl -lpthread -ldl -o "${FUNC}_generator"
./"${FUNC}_generator" > /dev/null

# ---------------------------------------------------------------------------
# 2. Compile the standalone measurement wrapper
# ---------------------------------------------------------------------------
echo "[2/3] Compiling the measurement wrapper"
${CXX} ${CXXFLAGS} ${INCLUDES} "${FUNC}_wrapper.cpp" "${FUNC}.o" ${LIB_FLAGS} \
    -ltiramisu -lHalide -lpthread -ldl -lm -o "${FUNC}_run"

# ---------------------------------------------------------------------------
# 3. Measure (PolyBench methodology: one kernel call per fresh process)
# ---------------------------------------------------------------------------
echo "[3/3] Running ${RUNS} measurement(s)"
if [ "${TPB_DUMP}" = "1" ]; then
    TIMES=$(TIRALIB_DUMP_ARRAYS=1 ./"${FUNC}_run" "${RUNS}" 2> "${FUNC}.dump")
    echo "live-out arrays dumped to ${KERNEL_FOLDER}/${FUNC}.dump"
else
    TIMES=$(./"${FUNC}_run" "${RUNS}")
fi

echo "times (ms): ${TIMES}"
if command -v python3 > /dev/null; then
    echo "${TIMES}" | python3 -c '
import sys
times = [float(x) for x in sys.stdin.read().split()]
if len(times) >= 3:
    mid = sorted(times)[1:-1]
    mean = sum(mid) / len(mid)
    dev = max(abs(t - mean) for t in mid) / mean * 100 if mean else 0.0
    print(f"polybench-normalized time (drop min/max, mean of middle {len(mid)}): "
          f"{mean:.6f} ms (max deviation {dev:.2f}%)")
'
fi

if [ "${TPB_KEEP}" != "1" ]; then
    rm -f "${FUNC}_generator" "${FUNC}_run" "${FUNC}.o" "${FUNC}.o.h"
fi
