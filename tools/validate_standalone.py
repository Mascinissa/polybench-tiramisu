#!/usr/bin/env python3
"""Validate the standalone path against stock PolyBench/C.

For each benchmark (default: all 30, MINI dataset):

1. Builds and runs the *stock* PolyBench/C benchmark
   (``-DPOLYBENCH_DUMP_ARRAYS``): init -> kernel -> live-out dump.
2. Runs ``compile_and_run.sh <bench> <SIZE> 1`` with ``TPB_DUMP=1``:
   Tiramisu generator -> standalone wrapper -> sidecar init -> kernel ->
   live-out dump.
3. Compares the two dumps numerically (same comparison as
   validate_against_polybench_c.py, including the gramschmidt
   singular-pivot masking).

Requires a Tiramisu installation (TIRAMISU_ROOT / TPB_* env, see
compile_and_run.sh) and a PolyBench/C 4.2.1 checkout — but *not*
TiraLib: this validates the pure standalone path.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from polybench_extract import BENCHMARKS, DATASETS  # noqa: E402
from validate_against_polybench_c import compare  # noqa: E402


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
    parser.add_argument("--dataset", default="MINI")
    parser.add_argument("--benchmarks", nargs="*", default=None)
    parser.add_argument("--atol", type=float, default=0.011)
    parser.add_argument("--rtol", type=float, default=1e-4)
    args = parser.parse_args()

    benchmarks = args.benchmarks or sorted(BENCHMARKS)
    results = []

    with tempfile.TemporaryDirectory(prefix="tpb_standalone_") as tmp:
        tmpdir = Path(tmp)
        for bench in benchmarks:
            name = f"function_{bench}_{args.dataset}"
            try:
                # --- stock reference ---
                stock_c = args.polybench_c / (BENCHMARKS[bench] + ".c")
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

                # --- standalone path ---
                env = dict(os.environ)
                env["TPB_DUMP"] = "1"
                build = subprocess.run(
                    [str(args.repo_root / "compile_and_run.sh"),
                     bench, args.dataset, "1"],
                    env=env, capture_output=True, text=True, check=True,
                )
                dumps = list(args.repo_root.glob(f"**/{name}.dump"))
                assert len(dumps) == 1, (dumps, build.stdout[-500:])
                standalone_dump = dumps[0].read_text()
                dumps[0].unlink()

                err, note = compare(
                    bench, standalone_dump, ref.stderr, args.atol, args.rtol
                )
                results.append((name, err or f"OK{note}"))
            except subprocess.CalledProcessError as e:
                detail = ((e.stderr or "") + (e.stdout or ""))[-300:]
                results.append((name, f"build/run failed: {detail}"))
            except Exception as e:  # noqa: BLE001
                results.append((name, f"error: {type(e).__name__}: {e}"))
            print(f"{results[-1][0]:45s} {results[-1][1]}", flush=True)

    n_ok = sum(1 for _, r in results if r.startswith("OK"))
    print(f"\n{n_ok}/{len(results)} standalone runs match PolyBench/C")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
