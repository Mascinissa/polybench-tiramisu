#!/bin/bash

if [ $# -eq 0 ]; then
      echo "Usage: ./complile_and_run.sh <benchmark_name> <problem_size>"
      echo "Example: ./complile_and_run.sh atax MEDIUM"
      exit
fi

if [ -z ${TIRAMISU_ROOT} ]
then
      echo "Tiramisu path not defined. Please specify the path to the Tiramisu root folder"
      echo "export TIRAMISU_ROOT=<path/to/Tiramisu>"
      exit
fi

# Define data sizes, possible value: -DTIRAMISU_XLARGE, -DTIRAMISU_LARGE, -DTIRAMISU_MEDIUM, -DTIRAMISU_SMALL
PBSIZE=$2
if [[ $2 = "MINI" ]]; 
then
      DEFINED_SIZE="-DTIRAMISU_MINI"
elif [[ $2 = "SMALL" ]]
then
      DEFINED_SIZE="-DTIRAMISU_SMALL"
elif [[ $2 = "MEDIUM" ]]
then
      DEFINED_SIZE="-DTIRAMISU_MEDIUM"
elif [[ $2 = "LARGE" ]]
then
      DEFINED_SIZE="-DTIRAMISU_LARGE"
elif [[ $2 = "EXTRALARGE" ]]
then
      DEFINED_SIZE="-DTIRAMISU_XLARGE"
else 
      echo "Unrecognized size"
      PBSIZE="MEDIUM"
      DEFINED_SIZE="-DTIRAMISU_MEDIUM"
fi

KERNEL=$1

if [ ${KERNEL} = "correlation" ]; then
      KERNEL_FOLDER="./datamining/correlation"
elif [ ${KERNEL} = "covariance" ]; then
      KERNEL_FOLDER="./datamining/covariance"
elif [ ${KERNEL} = "2mm" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/2mm"
elif [ ${KERNEL} = "3mm" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/3mm"
elif [ ${KERNEL} = "atax" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/atax"
elif [ ${KERNEL} = "bicg" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/bicg"
elif [ ${KERNEL} = "doitgen" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/doitgen"
elif [ ${KERNEL} = "mvt" ]; then
      KERNEL_FOLDER="./linear-algebra/kernels/mvt"
elif [ ${KERNEL} = "gemm" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/gemm"
elif [ ${KERNEL} = "gemver" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/gemver"
elif [ ${KERNEL} = "gesummv" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/gesummv"
elif [ ${KERNEL} = "symm" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/symm"
      echo "The symm benchmark has a known code generation issue (358aeab). Edit the compilation script to run it anyway."
      exit
elif [ ${KERNEL} = "syr2k" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/syr2k"
elif [ ${KERNEL} = "syrk" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/syrk"
elif [ ${KERNEL} = "trmm" ]; then
      KERNEL_FOLDER="./linear-algebra/blas/trmm"
elif [ ${KERNEL} = "cholesky" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/cholesky"
elif [ ${KERNEL} = "durbin" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/durbin"
elif [ ${KERNEL} = "gramschmidt" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/gramschmidt"
elif [ ${KERNEL} = "lu" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/lu"
elif [ ${KERNEL} = "ludcmp" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/ludcmp"
elif [ ${KERNEL} = "trisolv" ]; then
      KERNEL_FOLDER="./linear-algebra/solvers/trisolv"
elif [ ${KERNEL} = "deriche" ]; then
      KERNEL_FOLDER="./medley/deriche"
elif [ ${KERNEL} = "floyd_warshall" ]; then
      KERNEL_FOLDER="./medley/floyd_warshall"
elif [ ${KERNEL} = "nussinov" ]; then
      KERNEL_FOLDER="./medley/nussinov"
elif [ ${KERNEL} = "adi" ]; then
      KERNEL_FOLDER="./stencils/adi"
elif [ ${KERNEL} = "fdtd_2d" ]; then
      KERNEL_FOLDER="./stencils/fdtd_2d"
elif [ ${KERNEL} = "heat_3d" ]; then
      KERNEL_FOLDER="./stencils/heat_3d"
elif [ ${KERNEL} = "jacobi_1d" ]; then
      KERNEL_FOLDER="./stencils/jacobi_1d"
elif [ ${KERNEL} = "jacobi_2d" ]; then
      KERNEL_FOLDER="./stencils/jacobi_2d"
elif [ ${KERNEL} = "seidel_2d" ]; then
      KERNEL_FOLDER="./stencils/seidel_2d"
else 
      echo "Unrecognized benchmark"
      exit
fi

CORES=4
EXTRA_LIBRARIES="-ldl"

# Paths to Tiramisu 3rd party libraries
ISL_INCLUDE_DIRECTORY=${TIRAMISU_ROOT}/3rdParty/isl/build/include/
ISL_LIB_DIRECTORY=${TIRAMISU_ROOT}/3rdParty/isl/build/lib/
HALIDE_SOURCE_DIRECTORY=${TIRAMISU_ROOT}/3rdParty/Halide
HALIDE_LIB_DIRECTORY=${TIRAMISU_ROOT}/3rdParty/Halide/lib

CXXFLAGS="-std=c++11 -O3 -fno-rtti -mavx2"

CXX=g++

# Compile options
# - Make ${CXX} dump generated assembly
#   CXXFLAGS: -g -Wa,-alh
# - Get info about ${CXX} vectorization
#   CXXFLAGS -fopt-info-vec
# - Pass options to the llvm compiler
#   HL_LLVM_ARGS="-help" 
# - Set thread number for Halide
#   HL_NUM_THREADS=32
# - Execution env variables
#   OMP_NUM_THREADS=48
#   to set the number of threads to use by OpenMP.
# - Command to run Vtune
#   source /data/scratch/yunming/intel_parallel_studio_cluster/parallel_studio_xe_2017/install/vtune_amplifier_xe/amplxe-vars.sh
#   amplxe-cl -collect hpc-performance -result-dir vtune_results -quiet ./binary
#   Guide: https://software.intel.com/en-us/vtune-amplifier-help-amplxe-cl-command-syntax

INCLUDES="-I${MKL_PREFIX}/include/ -I${TIRAMISU_ROOT}/include/ -I${HALIDE_SOURCE_DIRECTORY}/include/ -I${ISL_INCLUDE_DIRECTORY} -I${KERNEL_FOLDER}/ -I${PWD}/utilities/"
LIBRARIES="-ltiramisu ${MKL_FLAGS} -lHalide -lisl -lz -lpthread ${EXTRA_LIBRARIES}"
LIBRARIES_DIR="-L${MKL_PREFIX}/lib/${MKL_LIB_PATH_SUFFIX} -L${HALIDE_LIB_DIRECTORY}/ -L${ISL_LIB_DIRECTORY}/ -L${TIRAMISU_ROOT}/build/"

echo "Compiling ${KERNEL} ${PBSIZE} "

cd ${KERNEL_FOLDER}

rm -rf ${KERNEL}_generator ${KERNEL}_wrapper generated_${KERNEL}.o generated_${KERNEL}_halide.o

# Generate code from Tiramisu
${CXX} ${LANKA_OPTIONS} $CXXFLAGS ${INCLUDES} ${DEFINED_SIZE} ${KERNEL}_generator.cpp ${LIBRARIES_DIR} ${LIBRARIES}                       -o ${KERNEL}_generator
echo "Running ${KERNEL} ${PBSIZE} generator (Tiramisu)"

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HALIDE_LIB_DIRECTORY}:${ISL_LIB_DIRECTORY}:${TIRAMISU_ROOT}/build/:${MKL_PREFIX}/lib/${MKL_LIB_PATH_SUFFIX} DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${HALIDE_LIB_DIRECTORY}:${TIRAMISU_ROOT}/build/:${MKL_PREFIX}/lib/${MKL_LIB_PATH_SUFFIX} ./${KERNEL}_generator

if [ $? -ne 0 ]; then
	exit
fi

# echo "Compiling ${KERNEL} wrapper"
${CXX} ${LANKA_OPTIONS} $CXXFLAGS ${INCLUDES} ${DEFINED_SIZE} ${KERNEL}_wrapper.cpp   ${LIBRARIES_DIR} ${LIBRARIES} generated_${KERNEL}.o ${LIBRARIES} -o ${KERNEL}_wrapper
echo "Running ${KERNEL} ${PBSIZE} wrapper"
# To enable profiling:
## Perf:
#PROFILING_COMMAND="perf stat -e cycles,instructions,cache-misses,L1-icache-load-misses,LLC-load-misses,dTLB-load-misses,cpu-migrations,context-switches,bus-cycles,cache-references,minor-faults"
## Vtune:
#VTUNE_METRIC=hpc-performance
#VTUNE_METRIC=memory-access
#PROFILING_COMMAND="amplxe-cl -collect ${VTUNE_METRIC} -result-dir vtune_results -quiet"
#rm -rf vtune_results

RUN_REF=1 RUN_TIRAMISU=1 HL_NUM_THREADS=$CORES LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HALIDE_LIB_DIRECTORY}:${ISL_LIB_DIRECTORY}:${TIRAMISU_ROOT}/build/:${MKL_PREFIX}/lib/${MKL_LIB_PATH_SUFFIX} DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${HALIDE_LIB_DIRECTORY}:${TIRAMISU_ROOT}/build/:${MKL_PREFIX}/lib/${MKL_LIB_PATH_SUFFIX} ${PROFILING_COMMAND} ./${KERNEL}_wrapper

rm -rf ${KERNEL}_generator ${KERNEL}_wrapper generated_${KERNEL}.o generated_${KERNEL}.o.h

cd -
