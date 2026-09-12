#!/bin/bash
# Fetch PolyBench/C 4.2.1 (needed only for the tools: sidecar
# regeneration and the parity/end-to-end validators — NOT for building
# or running the benchmarks, which use the vendored utilities/polybench.{c,h}).
#
# Usage: ./tools/get_polybench_c.sh [dest_dir]   (default: ../PolyBenchC-4.2.1)
set -e

DEST=${1:-"$(cd "$(dirname "$0")/.." && pwd)/../PolyBenchC-4.2.1"}
if [ -f "${DEST}/utilities/polybench.c" ]; then
    echo "PolyBench/C already present at ${DEST}"
    exit 0
fi

TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT

URL="https://downloads.sourceforge.net/project/polybench/polybench-c-4.2.1-beta.tar.gz"
echo "Downloading PolyBench/C 4.2.1 from ${URL} ..."
if ! curl -fL --retry 3 -o "${TMP}/polybench.tar.gz" "${URL}"; then
    echo "Download failed. Fetch PolyBench/C 4.2.1 manually (polybench.sourceforge.net)"
    echo "and point the tools at it with --polybench-c <path>."
    exit 1
fi

mkdir -p "${DEST}"
tar -xzf "${TMP}/polybench.tar.gz" -C "${TMP}"
EXTRACTED=$(find "${TMP}" -maxdepth 1 -mindepth 1 -type d | head -n 1)
cp -r "${EXTRACTED}"/. "${DEST}"/

if [ ! -f "${DEST}/utilities/polybench.c" ]; then
    echo "Unexpected archive layout under ${DEST}; check manually."
    exit 1
fi
echo "PolyBench/C 4.2.1 available at ${DEST}"
