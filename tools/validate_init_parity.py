#!/usr/bin/env python3
"""Validate the generated init sidecars against stock PolyBench/C.

For each ``function_<bench>_<SIZE>_init.h`` sidecar this script builds and
runs two small programs and compares their array dumps numerically:

1. **Sidecar program** (C++): PolyBench's polybench.h/polybench.c inlined
   (exactly as TiraLib's PolybenchHarness does), then the sidecar, then a
   main() that allocates the flat buffers, calls
   ``tiralib_polybench_init_arrays`` and ``tiralib_polybench_dump_arrays``.

2. **Reference program** (C): includes the *stock* PolyBench benchmark
   source (with its ``main`` renamed away), declares the arrays with
   PolyBench's own ``POLYBENCH_*_ARRAY_DECL`` macros, calls the stock
   ``init_array`` and ``print_array``. Arrays that ``init_array`` leaves
   untouched are zeroed (matching the sidecar's deterministic zero-fill).

A sidecar passes if both dumps contain the same sequence of numbers
(within a small tolerance that absorbs int/float/double formatting
differences). This validates the verbatim body extraction, the
flat-buffer casts, the baked dataset sizes, the buffer<->array mapping,
and the zero-fill decisions — everything the generator could get wrong.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from polybench_extract import (  # noqa: E402
    DATASETS,
    StockBenchmark,
    array_written_in_body,
    parse_stock_benchmark,
    parse_tiramisu_variant,
    stock_array_name,
)

NUM_RE = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


def extract_numbers(dump: str) -> list[float]:
    return [float(x) for x in NUM_RE.findall(dump)]


def build_sidecar_program(
    polybench_root: Path, sidecar_path: Path, tv
) -> str:
    polybench_h = (polybench_root / "utilities/polybench.h").read_text()
    polybench_c = (polybench_root / "utilities/polybench.c").read_text()
    polybench_c = re.sub(r'#\s*include\s+"polybench\.h"', "", polybench_c)

    decls, args, frees = [], [], []
    for buf, dims, ctype in zip(
        tv.buffer_names, tv.buffer_dims, tv.buffer_ctypes
    ):
        size = "*".join(str(d) for d in dims)
        decls.append(
            f"    {ctype} *{buf} = ({ctype}*) polybench_alloc_data"
            f"({size}, sizeof({ctype}));"
        )
        args.append(buf)
        frees.append(f"    free((void*){buf});")

    nl = "\n"
    return f"""\
#define POLYBENCH_TIME 1
{polybench_h}
{polybench_c}
#include "{sidecar_path}"

int main() {{
{nl.join(decls)}
    tiralib_polybench_init_arrays({", ".join(args)});
    tiralib_polybench_dump_arrays({", ".join(args)});
{nl.join(frees)}
    return 0;
}}
"""


def build_reference_program(stock: StockBenchmark, dataset: str) -> str:
    """dataset: PolyBench dataset name (e.g. EXTRALARGE)."""
    # Union of parameters over init_array and print_array, keeping order.
    int_params: list[str] = []
    scalar_params: list[str] = []
    arrays: dict[str, object] = {}
    for func in (stock.init, stock.print_):
        for p in func.int_params:
            if p not in int_params:
                int_params.append(p)
        for s in func.scalar_params:
            if s not in scalar_params:
                scalar_params.append(s)
        for arr in func.arrays:
            arrays.setdefault(arr.name, arr)

    lines = []
    for p in int_params:
        lines.append(f"    int {p} = {p.upper()};")
    for s in scalar_params:
        lines.append(f"    DATA_TYPE {s};")
    for arr in arrays.values():
        k = len(arr.macro_dims)
        macro_args = ", ".join(arr.macro_dims)
        lower_args = ", ".join(arr.lower_dims)
        lines.append(
            f"    POLYBENCH_{k}D_ARRAY_DECL({arr.name}, {arr.elem_type}, "
            f"{macro_args}, {lower_args});"
        )
    # Zero arrays the stock init does not write (the sidecar zero-fills
    # its corresponding buffers).
    for arr in arrays.values():
        if not array_written_in_body(stock.init.body, arr.name):
            lines.append(
                f"    memset((void*){arr.name}, 0, sizeof(*{arr.name}));"
            )

    def call_args(func) -> str:
        parts = list(func.int_params)
        parts += [f"&{s}" for s in func.scalar_params]
        parts += [f"POLYBENCH_ARRAY({a.name})" for a in func.arrays]
        return ", ".join(parts)

    lines.append(f"    init_array({call_args(stock.init)});")
    lines.append(f"    print_array({call_args(stock.print_)});")

    nl = "\n"
    return f"""\
#define {dataset}_DATASET 1
#define POLYBENCH_DUMP_ARRAYS 1
#define main tiralib_disabled_main
#include "{stock.c_path}"
#undef main
#include <string.h>

int main() {{
{nl.join(lines)}
    return 0;
}}
"""


def run_and_capture_stderr(cmd: list[str], cwd: Path) -> str:
    res = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, check=True
    )
    return res.stderr


def compare(sidecar_nums: list[float], ref_nums: list[float]) -> str | None:
    if len(sidecar_nums) != len(ref_nums):
        return (
            f"dump length mismatch: sidecar={len(sidecar_nums)} "
            f"ref={len(ref_nums)}"
        )
    for i, (a, b) in enumerate(zip(sidecar_nums, ref_nums)):
        if abs(a - b) > 0.02 and abs(a - b) > 1e-6 * max(abs(a), abs(b)):
            return f"value mismatch at #{i}: sidecar={a} ref={b}"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--polybench-c",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent
        / "PolyBenchC-4.2.1",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=["MINI", "SMALL", "MEDIUM"],
        help="Tiramisu dataset suffixes to run the parity check for "
        "(default: MINI SMALL MEDIUM; larger dumps get slow)",
    )
    parser.add_argument(
        "--compile-only-datasets",
        nargs="*",
        default=["LARGE", "XLARGE"],
        help="Datasets for which the sidecar program is only compiled",
    )
    parser.add_argument("--benchmarks", nargs="*", default=None)
    parser.add_argument("--keep-going", action="store_true")
    args = parser.parse_args()

    generators = sorted(
        p for p in args.repo_root.glob("**/function_*_generator.cpp")
        if ".git" not in p.parts
    )
    stock_cache: dict[str, StockBenchmark] = {}
    passed, compiled, failures = 0, 0, []

    with tempfile.TemporaryDirectory(prefix="tiralib_parity_") as tmp:
        tmpdir = Path(tmp)
        for gen_path in generators:
            tv = parse_tiramisu_variant(gen_path)
            if args.benchmarks and tv.benchmark not in args.benchmarks:
                continue
            run_parity = tv.dataset in args.datasets
            compile_only = tv.dataset in args.compile_only_datasets
            if not (run_parity or compile_only):
                continue
            sidecar_path = gen_path.parent / f"{tv.function_name}_init.h"
            if not sidecar_path.exists():
                failures.append((tv.function_name, "sidecar missing"))
                continue
            if tv.benchmark not in stock_cache:
                stock_cache[tv.benchmark] = parse_stock_benchmark(
                    args.polybench_c, tv.benchmark
                )
            stock = stock_cache[tv.benchmark]

            try:
                # --- sidecar program (C++, PolybenchHarness-style TU) ---
                sc_src = tmpdir / f"{tv.function_name}_sidecar.cpp"
                sc_bin = tmpdir / f"{tv.function_name}_sidecar"
                sc_src.write_text(
                    build_sidecar_program(
                        args.polybench_c, sidecar_path, tv
                    )
                )
                subprocess.run(
                    ["g++", "-std=c++17", "-O1", "-o", str(sc_bin),
                     str(sc_src), "-lm"],
                    check=True, capture_output=True, text=True,
                )
                compiled += 1
                if not run_parity:
                    continue
                sidecar_dump = run_and_capture_stderr([str(sc_bin)], tmpdir)

                # --- reference program (C, stock PolyBench) ---
                ref_src = tmpdir / f"{tv.function_name}_ref.c"
                ref_bin = tmpdir / f"{tv.function_name}_ref"
                ref_src.write_text(
                    build_reference_program(stock, DATASETS[tv.dataset])
                )
                subprocess.run(
                    ["gcc", "-O1",
                     "-I", str(args.polybench_c / "utilities"),
                     "-I", str(stock.c_path.parent),
                     str(args.polybench_c / "utilities/polybench.c"),
                     str(ref_src), "-o", str(ref_bin), "-lm"],
                    check=True, capture_output=True, text=True,
                )
                ref_dump = run_and_capture_stderr([str(ref_bin)], tmpdir)

                err = compare(
                    extract_numbers(sidecar_dump), extract_numbers(ref_dump)
                )
                if err:
                    failures.append((tv.function_name, err))
                    if not args.keep_going:
                        break
                else:
                    passed += 1
            except subprocess.CalledProcessError as e:
                failures.append(
                    (tv.function_name,
                     f"build/run failed: {e.stderr[:500] if e.stderr else e}")
                )
                if not args.keep_going:
                    break

    print(f"parity PASS: {passed}, compile-only OK: {compiled - passed}")
    if failures:
        print(f"{len(failures)} FAILURES:", file=sys.stderr)
        for name, err in failures:
            print(f"  {name}: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
