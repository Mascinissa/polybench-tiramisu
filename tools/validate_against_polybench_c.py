#!/usr/bin/env python3
"""End-to-end validation of the PolyBench harness against PolyBench/C.

For each benchmark (default: all 30, MINI dataset):

1. Builds and runs the *stock* PolyBench/C benchmark
   (``-DPOLYBENCH_DUMP_ARRAYS``): init -> kernel -> live-out dump.
2. Loads the Tiramisu generator through TiraLib (the PolyBench init
   sidecar is auto-detected, so the PolybenchHarness is used), executes
   an empty schedule, then re-runs the produced measurement wrapper with
   ``TIRALIB_DUMP_ARRAYS=1``: sidecar init -> tiramisu kernel ->
   live-out dump.
3. Compares the two dumps numerically.

A match means the whole chain is faithful: PolyBench initialization,
buffer layout/mapping, the Tiramisu reimplementation's semantics, and
the dump of the same live-out arrays.

Run inside the tiramisu conda environment.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from polybench_extract import (  # noqa: E402
    BENCHMARKS,
    DATASETS,
    parse_tiramisu_variant,
)

def extract_sections(dump: str) -> dict[str, list[float]]:
    """Parse a PolyBench POLYBENCH_DUMP_ARRAYS dump into per-array value
    lists. Tokens like nan/inf parse to the corresponding floats so both
    dumps stay aligned position-by-position."""
    sections: dict[str, list[float]] = {}
    for m in re.finditer(
        r"begin dump: (\w+)(.*?)end   dump", dump, re.DOTALL
    ):
        vals = []
        for tok in m.group(2).split():
            try:
                vals.append(float(tok))
            except ValueError:
                pass
        sections[m.group(1)] = vals
    return sections


def compare_values(a, b, atol, rtol, mask=None, label=""):
    if len(a) != len(b):
        return f"{label}: length mismatch: tiramisu={len(a)} stock={len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        if mask is not None and mask(i):
            continue
        if math.isnan(x) != math.isnan(y):
            return f"{label}[{i}]: nan mismatch: tiramisu={x} stock={y}"
        if math.isnan(x):
            continue
        d = abs(x - y)
        if d > atol and d > rtol * max(abs(x), abs(y)):
            return f"{label}[{i}]: tiramisu={x} stock={y} (absdiff={d})"
    return None


def compare(bench, tiramisu_dump, stock_dump, atol, rtol):
    """Compare two PolyBench dumps. Returns (error, note)."""
    t = extract_sections(tiramisu_dump)
    s = extract_sections(stock_dump)
    if list(t) != list(s):
        return f"dumped arrays differ: tiramisu={list(t)} stock={list(s)}", ""

    note = ""
    masks = {name: None for name in s}
    if bench == "gramschmidt" and "R" in s:
        # PolyBench's gramschmidt input is rank-deficient (N > M and the
        # modular init), so classical Gram-Schmidt hits a numerically
        # zero pivot R[k0][k0]; beyond it, Q[:,k>=k0] and R[k>=k0][:] are
        # normalized cancellation residue (~1e-15) whose value is
        # FP-contraction dependent — the stock binary itself differs
        # there when compiled with -mfma -ffp-contract=fast. Compare only
        # the well-defined prefix.
        n = math.isqrt(len(s["R"]))
        k0 = next(
            (k for k in range(n) if abs(s["R"][k * n + k]) <= atol), n
        )
        if k0 < n:
            note = f" (up to singular pivot k={k0}/{n}; tail is FP-undefined)"
            masks["R"] = lambda i, n=n, k0=k0: i // n >= k0
            if "Q" in s:
                masks["Q"] = lambda i, n=n, k0=k0: i % n >= k0

    for name in s:
        err = compare_values(
            t[name], s[name], atol, rtol, mask=masks.get(name), label=name
        )
        if err:
            return err, note
    return None, note


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--polybench-c", type=Path,
        default=Path(__file__).resolve().parent.parent.parent
        / "PolyBenchC-4.2.1",
    )
    parser.add_argument(
        "--repo-root", type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument(
        "--tiralib-root", type=Path,
        default=Path(__file__).resolve().parent.parent.parent / "TiraLib",
    )
    parser.add_argument("--dataset", default="MINI",
                        help="Tiramisu dataset suffix (default MINI)")
    parser.add_argument("--benchmarks", nargs="*", default=None)
    parser.add_argument("--atol", type=float, default=0.011,
                        help="Absolute tolerance (dump print quantum)")
    parser.add_argument("--rtol", type=float, default=1e-4)
    args = parser.parse_args()

    sys.path.insert(0, str(args.tiralib_root))
    os.chdir(args.tiralib_root)  # so BaseConfig.init finds config.yaml

    from tiralib.config import BaseConfig
    from tiralib.tiramisu.harness import PolybenchHarness
    from tiralib.tiramisu.schedule import Schedule
    from tiralib.tiramisu.tiramisu_program import TiramisuProgram

    import logging
    BaseConfig.init(logging_level=logging.WARNING)
    config = BaseConfig.base_config
    workspace = Path(config.workspace)

    # Environment for running the wrapper manually (mirrors
    # CompilingService.get_env_vars).
    wrapper_env = dict(os.environ)
    libs = ":".join(config.dependencies.libs)
    wrapper_env["LD_LIBRARY_PATH"] = (
        libs + ":" + wrapper_env.get("LD_LIBRARY_PATH", "")
    )

    benchmarks = args.benchmarks or sorted(BENCHMARKS)
    results: list[tuple[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="tiralib_e2e_") as tmp:
        tmpdir = Path(tmp)
        generator_index = {
            g.stem[: -len("_generator")]: g
            for g in args.repo_root.glob("**/function_*_generator.cpp")
            if ".git" not in g.parts
        }
        for bench in benchmarks:
            name = f"function_{bench}_{args.dataset}"
            gen_path = generator_index.get(name)
            if gen_path is None:
                results.append((name, "generator not found in repo tree"))
                print(f"{results[-1][0]:45s} {results[-1][1]}", flush=True)
                continue
            try:
                tv = parse_tiramisu_variant(gen_path)

                # --- stock reference: init -> kernel -> dump ---
                stock_rel = BENCHMARKS[bench]
                stock_c = args.polybench_c / (stock_rel + ".c")
                ref_bin = tmpdir / f"{bench}_ref"
                subprocess.run(
                    ["gcc", "-O2",
                     "-I", str(args.polybench_c / "utilities"),
                     "-I", str(stock_c.parent),
                     f"-D{DATASETS[args.dataset]}_DATASET",
                     "-DPOLYBENCH_DUMP_ARRAYS",
                     str(args.polybench_c / "utilities/polybench.c"),
                     str(stock_c), "-o", str(ref_bin), "-lm"],
                    check=True, capture_output=True, text=True,
                )
                ref = subprocess.run(
                    [str(ref_bin)], capture_output=True, text=True,
                    check=True, cwd=tmpdir,
                )

                # --- TiraLib with PolybenchHarness ---
                program = TiramisuProgram.from_file(
                    str(gen_path), load_annotations=True, load_tree=True
                )
                assert isinstance(program.harness, PolybenchHarness), (
                    f"{name}: PolybenchHarness not auto-detected"
                )
                schedule = Schedule(program)
                times = schedule.execute(min_runs=1, delete_files=False)
                assert len(times) == 1 and times[0] > 0, times

                run_env = dict(wrapper_env)
                run_env["MIN_RUNS"] = "1"
                run_env["TIRALIB_DUMP_ARRAYS"] = "1"
                wrapper = subprocess.run(
                    [f"./{program.temp_files_identifier}_wrapper"],
                    cwd=workspace, env=run_env,
                    capture_output=True, text=True, check=True,
                )
                err, note = compare(
                    bench, wrapper.stderr, ref.stderr, args.atol, args.rtol
                )
                results.append((name, err or f"OK{note}"))
            except AssertionError as e:
                results.append((name, f"ASSERT: {e}"))
            except subprocess.CalledProcessError as e:
                detail = (e.stderr or "")[-300:]
                results.append((name, f"build/run failed: {detail}"))
            except Exception as e:  # noqa: BLE001
                results.append((name, f"error: {type(e).__name__}: {e}"))
            finally:
                # clean the TiraLib workspace files of this program
                try:
                    for f in workspace.glob(
                        f"{program.temp_files_identifier}*"
                    ):
                        f.unlink()
                except Exception:  # noqa: BLE001
                    pass

            print(f"{results[-1][0]:45s} {results[-1][1]}", flush=True)

    n_ok = sum(1 for _, r in results if r.startswith("OK"))
    print(f"\n{n_ok}/{len(results)} benchmarks match PolyBench/C")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
