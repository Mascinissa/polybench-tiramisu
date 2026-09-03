#!/usr/bin/env python3
"""Generate the per-benchmark init sidecars and standalone wrappers.

For every ``function_<bench>_<SIZE>`` variant in this repository (found
recursively), this script generates ``function_<bench>_<SIZE>_init.h``
(the init sidecar, auto-detected by TiraLib's PolybenchHarness) and
``function_<bench>_<SIZE>_wrapper.cpp``/``_wrapper.h`` (the standalone,
TiraLib-free measurement wrapper built by compile_and_run.sh) next to
the generator file. The sidecar contains the benchmark's *verbatim* ``init_array`` and
``print_array`` bodies extracted from PolyBench/C 4.2.1, wrapped in two
functions operating on the flat buffers that the Tiramisu function takes
(in ``tiramisu::codegen`` order)::

    static void tiralib_polybench_init_arrays(<typed buffer ptrs>);
    static void tiralib_polybench_dump_arrays(<typed buffer ptrs>);

TiraLib's ``PolybenchHarness`` auto-detects these sidecars and inlines
them into its measurement wrapper (after PolyBench's own polybench.h/.c),
so PolyBench macros like ``POLYBENCH_2D_ARRAY_DECL`` used inside some
init bodies (cholesky, lu, ...) resolve to PolyBench's own definitions.

Exposed buffers that PolyBench's ``init_array`` does not initialize
(kernel-written live-out arrays such as atax's ``y``) are zero-filled,
deterministically.

The script hard-fails on any mapping or dimension inconsistency between
the Tiramisu buffers and the stock PolyBench arrays.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from polybench_extract import (  # noqa: E402
    DATA_TYPE_INFO,
    DATASETS,
    FuncInfo,
    StockBenchmark,
    TiramisuVariant,
    array_written_in_body,
    parse_stock_benchmark,
    parse_tiramisu_variant,
    stock_array_name,
)


def _buffer_index(tv: TiramisuVariant, sname: str) -> int | None:
    """Index of the buffer mapping to stock array `sname`, or None."""
    for i, buf in enumerate(tv.buffer_names):
        if stock_array_name(tv.benchmark, buf) == sname:
            return i
    return None


def _check_dims(
    stock: StockBenchmark,
    tv: TiramisuVariant,
    func: FuncInfo,
    func_name: str,
    ds_macros: dict[str, int],
) -> None:
    for arr in func.arrays:
        idx = _buffer_index(tv, arr.name)
        if idx is None:
            raise ValueError(
                f"{tv.function_name}: stock {func_name} array {arr.name!r} "
                f"has no matching Tiramisu buffer "
                f"(buffers: {tv.buffer_names})"
            )
        macro_vals = [ds_macros[d] for d in arr.macro_dims]
        if macro_vals != tv.buffer_dims[idx]:
            raise ValueError(
                f"{tv.function_name}: dimension mismatch for {arr.name}: "
                f"stock {arr.macro_dims}={macro_vals} vs tiramisu buffer "
                f"{tv.buffer_names[idx]}={tv.buffer_dims[idx]}"
            )


def _func_prelude(
    stock: StockBenchmark,
    tv: TiramisuVariant,
    func: FuncInfo,
    ds_macros: dict[str, int],
) -> tuple[list[str], set[int]]:
    """Declaration lines adapting the verbatim body to flat buffers.

    Returns (lines, used_buffer_indices).
    """
    lines: list[str] = []
    used: set[int] = set()

    for p in func.int_params:
        macro = p.upper()
        if macro not in ds_macros:
            raise ValueError(
                f"{tv.function_name}: size parameter {p!r} has no dataset "
                f"macro {macro} (have: {sorted(ds_macros)})"
            )
        lines.append(f"  const int {p} = {macro}; (void){p};")

    for s in func.scalar_params:
        lines.append(
            f"  DATA_TYPE _tl_scalar_{s}; "
            f"DATA_TYPE *{s} = &_tl_scalar_{s}; (void){s};"
        )

    for arr in func.arrays:
        idx = _buffer_index(tv, arr.name)
        assert idx is not None  # _check_dims ran before
        used.add(idx)
        buf = tv.buffer_names[idx]
        etype = arr.elem_type  # DATA_TYPE or a typedef (e.g. base)
        if len(arr.macro_dims) == 1:
            lines.append(f"  {etype} *{arr.name} = _tl_{buf}; (void){arr.name};")
        else:
            inner = "".join(f"[{d}]" for d in arr.macro_dims[1:])
            lines.append(
                f"  {etype} (*{arr.name}){inner} = "
                f"({etype} (*){inner}) _tl_{buf}; (void){arr.name};"
            )

    return lines, used


def build_sidecar(stock: StockBenchmark, tv: TiramisuVariant) -> str:
    ds_macros = stock.datasets[DATASETS[tv.dataset]]
    type_info = DATA_TYPE_INFO[stock.data_type]

    # --- consistency checks -------------------------------------------------
    _check_dims(stock, tv, stock.init, "init_array", ds_macros)
    _check_dims(stock, tv, stock.print_, "print_array", ds_macros)

    init_covered: set[int] = set()
    for arr in stock.init.arrays:
        idx = _buffer_index(tv, arr.name)
        assert idx is not None
        if arr.elem_type == "DATA_TYPE":
            if type_info["ctype"] != tv.buffer_ctypes[idx]:
                raise ValueError(
                    f"{tv.function_name}: type mismatch for {arr.name}: "
                    f"stock {type_info['ctype']} vs tiramisu "
                    f"{tv.buffer_ctypes[idx]}"
                )
        if array_written_in_body(stock.init.body, arr.name):
            init_covered.add(idx)

    # --- typedefs used in signatures (e.g. nussinov's `base`) --------------
    typedef_lines: list[str] = []
    for arr in stock.init.arrays + stock.print_.arrays:
        if arr.elem_type != "DATA_TYPE":
            if arr.elem_type not in stock.typedefs:
                raise ValueError(
                    f"{tv.function_name}: unknown array element type "
                    f"{arr.elem_type!r}"
                )
            idx = _buffer_index(tv, arr.name)
            assert idx is not None
            line = (
                f"typedef {tv.buffer_ctypes[idx]} {arr.elem_type}; "
                f"/* stock: typedef {stock.typedefs[arr.elem_type]} "
                f"{arr.elem_type}; value range preserved */"
            )
            if line not in typedef_lines:
                typedef_lines.append(line)

    # --- function signatures ------------------------------------------------
    params = ", ".join(
        f"{ctype} *_tl_{buf}"
        for buf, ctype in zip(tv.buffer_names, tv.buffer_ctypes)
    )

    # --- init function ------------------------------------------------------
    init_lines, init_used = _func_prelude(stock, tv, stock.init, ds_macros)
    zero_fill_lines: list[str] = []
    for i, buf in enumerate(tv.buffer_names):
        if i not in init_covered:
            size = "*".join(str(d) for d in tv.buffer_dims[i])
            sname = stock_array_name(tv.benchmark, buf)
            zero_fill_lines.append(
                f"  memset(_tl_{buf}, 0, ({size}) * "
                f"sizeof({tv.buffer_ctypes[i]})); "
                f"/* {sname}: not initialized by PolyBench init_array */"
            )
            init_used.add(i)
    for i, buf in enumerate(tv.buffer_names):
        if i not in init_used:
            init_lines.append(f"  (void)_tl_{buf};")

    # --- dump function ------------------------------------------------------
    dump_lines, dump_used = _func_prelude(stock, tv, stock.print_, ds_macros)
    for i, buf in enumerate(tv.buffer_names):
        if i not in dump_used:
            dump_lines.append(f"  (void)_tl_{buf};")

    guard = f"TIRALIB_POLYBENCH_INIT_{tv.function_name.upper()}_H"
    macro_defs = "\n".join(
        f"#define {m} {v}" for m, v in sorted(ds_macros.items())
    )
    macro_undefs = "\n".join(f"#undef {m}" for m in sorted(ds_macros))

    nl = "\n"
    return f"""\
/* Auto-generated by tools/generate_polybench_init_sidecars.py — DO NOT EDIT.
 *
 * PolyBench-compatible initialization + live-out dump for
 * {tv.function_name}, extracted verbatim from
 * PolyBenchC-4.2.1/{stock.c_path.parent.name}/{stock.c_path.name}
 * ({DATASETS[tv.dataset]}_DATASET).
 *
 * Contract: both functions take the Tiramisu function's I/O buffers as
 * flat pointers, in tiramisu::codegen order. This header must be
 * preceded by PolyBench's polybench.h/polybench.c in the same
 * translation unit (TiraLib's PolybenchHarness inlines them); some init
 * bodies use POLYBENCH_* allocation macros.
 */
#ifndef {guard}
#define {guard}

#include <math.h>
#include <stdio.h>
#include <string.h>

{macro_defs}

#define DATA_TYPE {type_info["ctype"]}
#define DATA_PRINTF_MODIFIER {type_info["printf"]}
#define SCALAR_VAL(x) {type_info["scalar_val"]}
#define SQRT_FUN(x) {type_info["sqrt"]}
#define EXP_FUN(x) {type_info["exp"]}
#define POW_FUN(x,y) {type_info["pow"]}
{nl.join(typedef_lines)}

static void tiralib_polybench_init_arrays({params})
{{
{nl.join(init_lines)}
{nl.join(zero_fill_lines)}
  /* begin verbatim PolyBench init_array body */
  {stock.init.body}
  /* end verbatim PolyBench init_array body */
}}

static void tiralib_polybench_dump_arrays({params})
{{
{nl.join(dump_lines)}
  /* begin verbatim PolyBench print_array body */
  {stock.print_.body}
  /* end verbatim PolyBench print_array body */
}}

{macro_undefs}
#undef DATA_TYPE
#undef DATA_PRINTF_MODIFIER
#undef SCALAR_VAL
#undef SQRT_FUN
#undef EXP_FUN
#undef POW_FUN

#endif /* {guard} */
"""



def build_standalone_wrapper_cpp(tv: TiramisuVariant) -> str:
    """Standalone (TiraLib-free) measurement wrapper for one variant.

    Same measurement semantics as the wrapper TiraLib's PolybenchHarness
    generates: PolyBench-aligned allocation, sidecar init, PolyBench
    timer around a single kernel call, milliseconds on stdout, optional
    live-out dump on stderr. The fresh-process-per-measurement driver
    lives in utilities/harness_main.inc.
    """
    decls, frees = [], []
    for buf, dims, ctype in zip(
        tv.buffer_names, tv.buffer_dims, tv.buffer_ctypes
    ):
        size = "*".join(str(d) for d in dims)
        halide_dims = ", ".join(str(d) for d in dims[::-1])
        decls.append(
            f"    {ctype} *c_{buf} = ({ctype}*) polybench_alloc_data"
            f"({size}, sizeof({ctype}));"
        )
        decls.append(
            f"    Halide::Buffer<{ctype}> {buf}(c_{buf}, {halide_dims});"
        )
        frees.append(f"    free((void*)c_{buf});")
    c_args = ", ".join(f"c_{b}" for b in tv.buffer_names)
    func_args = ", ".join(f"{b}.raw_buffer()" for b in tv.buffer_names)

    nl = "\n"
    return f"""\
/* Auto-generated by tools/generate_polybench_init_sidecars.py — DO NOT EDIT.
 * Standalone PolyBench-methodology measurement wrapper for
 * {tv.function_name}. Build via ./compile_and_run.sh; see
 * utilities/harness_main.inc for the execution protocol.
 */
#define POLYBENCH_TIME 1
#include "polybench.h"
#include "polybench.c"

#include "{tv.function_name}_init.h"

#include "Halide.h"
#include "{tv.function_name}_wrapper.h"
#include <cstdio>

static void tpb_single_measurement(int dump_arrays)
{{
{nl.join(decls)}

    tiralib_polybench_init_arrays({c_args});

    /* polybench_timer_start() == polybench_prepare_instruments()
     * (cache flush as configured) + timer start. */
    polybench_timer_start();
    {tv.function_name}({func_args});
    polybench_timer_stop();

    /* PolyBench prints seconds; the harness protocol is milliseconds. */
    std::printf("%0.6f ", (polybench_t_end - polybench_t_start) * 1000.0);
    std::fflush(stdout);

    if (dump_arrays) {{
        tiralib_polybench_dump_arrays({c_args});
        std::fflush(stderr);
    }}

{nl.join(frees)}
}}

#include "harness_main.inc"
"""


def build_standalone_wrapper_h(tv: TiramisuVariant) -> str:
    params = ", ".join(f"halide_buffer_t *{b}" for b in tv.buffer_names)
    return f"""\
/* Auto-generated by tools/generate_polybench_init_sidecars.py — DO NOT EDIT.
 * Declaration of the Tiramisu-generated kernel. Included both by the
 * generator (which requires it to exist) and by the wrapper.
 */
#include <tiramisu/utils.h>

#ifdef __cplusplus
extern "C" {{
#endif
int {tv.function_name}({params});
#ifdef __cplusplus
}}  // extern "C"
#endif
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--polybench-c",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent
        / "PolyBenchC-4.2.1",
        help="Path to the PolyBench/C 4.2.1 checkout",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Path to the polybench_tiramisu repository root",
    )
    parser.add_argument(
        "--benchmarks",
        nargs="*",
        default=None,
        help="Subset of benchmarks to generate (default: all)",
    )
    args = parser.parse_args()

    generators = sorted(
        p for p in args.repo_root.glob("**/function_*_generator.cpp")
        if ".git" not in p.parts
    )
    if not generators:
        print(f"No generators found under {args.repo_root}", file=sys.stderr)
        return 1

    stock_cache: dict[str, StockBenchmark] = {}
    written, failures = 0, []
    for gen_path in generators:
        try:
            tv = parse_tiramisu_variant(gen_path)
        except Exception as e:  # noqa: BLE001
            failures.append((gen_path.name, f"parse error: {e}"))
            continue
        if args.benchmarks and tv.benchmark not in args.benchmarks:
            continue
        try:
            if tv.benchmark not in stock_cache:
                stock_cache[tv.benchmark] = parse_stock_benchmark(
                    args.polybench_c, tv.benchmark
                )
            sidecar = build_sidecar(stock_cache[tv.benchmark], tv)
        except Exception as e:  # noqa: BLE001
            failures.append((tv.function_name, str(e)))
            continue
        out_path = gen_path.parent / f"{tv.function_name}_init.h"
        out_path.write_text(sidecar)
        (gen_path.parent / f"{tv.function_name}_wrapper.cpp").write_text(
            build_standalone_wrapper_cpp(tv)
        )
        (gen_path.parent / f"{tv.function_name}_wrapper.h").write_text(
            build_standalone_wrapper_h(tv)
        )
        written += 1

    print(f"Wrote {written} sidecars.")
    if failures:
        print(f"{len(failures)} FAILURES:", file=sys.stderr)
        for name, err in failures:
            print(f"  {name}: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
