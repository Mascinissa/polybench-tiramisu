#!/usr/bin/env python3
"""Package benchmarks into a flat archive for TiraLib-based pipelines.

Evaluation pipelines built on TiraLib (e.g. autotuning/genopt eval
scripts) consume benchmark sources from zip archives with the flat
layout::

    <function_name>/<function_name>_generator.cpp
    <function_name>/<function_name>_init.h

keyed by function name — independent of this repository's directory
taxonomy. This script walks the repo tree and produces such an archive.
Only the generator and the init sidecar are packaged (TiraLib generates
its own measurement wrapper from the sidecar); the standalone
``*_wrapper.cpp`` files are not included.

The archive also contains a ``MANIFEST.json`` (repo commit, generation
timestamp, benchmark/size list) so that evaluation results can be traced
back to the exact benchmark version.

By default the script refuses to package a variant whose init sidecar is
missing, so a PolyBench program can never silently fall back to
TiraLib's default (non-PolyBench) harness; pass ``--allow-missing-sidecars``
to override (e.g. for intentionally non-PolyBench programs).
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path


def repo_commit(repo_root: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True,
                        help="Output zip path")
    parser.add_argument(
        "--repo-root", type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument("--benchmarks", nargs="*", default=None,
                        help="Subset of benchmark names (default: all)")
    parser.add_argument("--sizes", nargs="*", default=None,
                        help="Subset of sizes, e.g. MINI MEDIUM (default: all)")
    parser.add_argument("--allow-missing-sidecars", action="store_true")
    args = parser.parse_args()

    generators = sorted(
        p for p in args.repo_root.glob("**/function_*_generator.cpp")
        if ".git" not in p.parts
    )

    packaged, errors = [], []
    with zipfile.ZipFile(args.out, "w", zipfile.ZIP_DEFLATED) as z:
        for gen in generators:
            name = gen.stem[: -len("_generator")]
            m = re.fullmatch(r"function_(\w+?)_([A-Z]+)", name)
            if not m:
                continue
            bench, size = m.group(1), m.group(2)
            if args.benchmarks and bench not in args.benchmarks:
                continue
            if args.sizes and size not in args.sizes:
                continue
            sidecar = gen.parent / f"{name}_init.h"
            if not sidecar.exists():
                if args.allow_missing_sidecars:
                    z.write(gen, f"{name}/{name}_generator.cpp")
                    packaged.append(name)
                    continue
                errors.append(f"{name}: init sidecar missing ({sidecar})")
                continue
            z.write(gen, f"{name}/{name}_generator.cpp")
            z.write(sidecar, f"{name}/{name}_init.h")
            packaged.append(name)

        manifest = {
            "repository": "polybench-tiramisu",
            "commit": repo_commit(args.repo_root),
            "generated_at": datetime.datetime.now(
                datetime.timezone.utc
            ).isoformat(),
            "variants": sorted(packaged),
        }
        z.writestr("MANIFEST.json", json.dumps(manifest, indent=2))

    if errors:
        print(f"{len(errors)} ERRORS (archive incomplete):", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        args.out.unlink(missing_ok=True)
        return 1

    print(f"Packaged {len(packaged)} variants into {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
