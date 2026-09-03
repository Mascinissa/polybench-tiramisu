# utilities/

- `polybench.h`, `polybench.c` — verbatim, unmodified copies from
  PolyBench/C 4.2.1 (see `NOTICE` at the repo root). Used by the
  standalone measurement wrappers; do not edit.
- `harness_main.inc` — shared driver for the generated standalone
  wrappers: fresh-process-per-measurement protocol, milliseconds on
  stdout, optional live-out dump on stderr.
