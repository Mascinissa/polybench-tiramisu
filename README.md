# polybench-tiramisu

PolyBench-Tiramisu is a reimplementation of [PolyBench/C 4.2.1](http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/) in the [Tiramisu](http://tiramisu-compiler.org/) polyhedral compiler. It is a benchmark suite of 30 numerical computations with static control flow, extracted from operations in various application domains (linear algebra, image processing, physics simulation, dynamic programming, statistics, ...).

**v2** is a major rework:

- All 30 benchmarks reimplemented as size-specialized Tiramisu programs (5 dataset sizes each, 150 variants), compatible with [TiraLib](https://github.com/Tiramisu-Compiler/TiraLib) for applying and evaluating schedules/transformations.
- **PolyBench-faithful measurement**: benchmarks are measured with PolyBench/C's own methodology — deterministic `init_array` initialization (extracted verbatim from PolyBench/C), page-aligned allocation, cache flush before the timer, `gettimeofday` timing, one kernel invocation per fresh process — using the vendored, unmodified `polybench.{c,h}`.
- **Correctness checking**: every variant can dump its live-out arrays in PolyBench's `POLYBENCH_DUMP_ARRAYS` format for numerical comparison against stock PolyBench/C. All 30 benchmarks are validated (see [Validation](#validation)).
- Works **standalone** (a Tiramisu installation and a C++ compiler; no TiraLib) and **with TiraLib** (its `PolybenchHarness` auto-detects these benchmarks — no glue code).

Looking for the original version? See the [`v1-legacy`](../../tree/v1-legacy) tag.

## Layout

```
linear-algebra/blas/gemm/
├── function_gemm_MINI_generator.cpp   # Tiramisu program (sizes baked in)
├── function_gemm_MINI_init.h          # init sidecar: verbatim PolyBench init_array/print_array (generated)
├── function_gemm_MINI_wrapper.cpp     # standalone measurement wrapper (generated)
├── function_gemm_MINI_wrapper.h
└── ... SMALL / MEDIUM / LARGE / XLARGE
utilities/          # vendored PolyBench/C polybench.{c,h} + shared measurement driver
tools/              # generation, validation, and packaging tools
```

## Standalone usage

Set `TIRAMISU_ROOT` to your Tiramisu installation (both the classic source layout and the `install/` layout are detected; `TPB_INCLUDES`/`TPB_LIBS` override detection — see the header of `compile_and_run.sh`).

```sh
export TIRAMISU_ROOT=<path/to/tiramisu>
./compile_and_run.sh gemm MEDIUM        # 5 measurements (PolyBench's count)
./compile_and_run.sh cholesky LARGE 10  # 10 measurements
```

Output: one time (milliseconds) per measurement, plus the PolyBench-normalized time (drop min and max, mean of the middle runs — `time_benchmark.sh`'s statistic):

```
times (ms): 4.640 4.634 4.617 4.710 4.727
polybench-normalized time (drop min/max, mean of middle 3): 4.661 ms (max deviation 1.04%)
```

Each measurement runs in a freshly exec'ed process and performs: aligned allocation → PolyBench `init_array` → cache flush + timer start → **one** kernel invocation → timer stop. This mirrors running a PolyBench/C binary N times.

To dump the live-out arrays (PolyBench `POLYBENCH_DUMP_ARRAYS` format, for correctness diffing):

```sh
TPB_DUMP=1 ./compile_and_run.sh gemm MINI 1
# -> linear-algebra/blas/gemm/function_gemm_MINI.dump
```

## Usage with TiraLib

The init sidecars are auto-detected by TiraLib, which switches to its `PolybenchHarness` (same measurement methodology, same wrapper semantics):

```python
from tiralib.config import BaseConfig
from tiralib.tiramisu import TiramisuProgram, Schedule, tiramisu_actions

BaseConfig.init()
program = TiramisuProgram.from_file(
    "linear-algebra/blas/gemm/function_gemm_MEDIUM_generator.cpp",
    load_annotations=True, load_tree=True,
)  # PolybenchHarness auto-detected via the *_init.h sidecar

schedule = Schedule(program)
schedule.add_optimizations([tiramisu_actions.Parallelization([("C_init", 0)])])
times_ms = schedule.execute(min_runs=5)   # 5 fresh-process PolyBench measurements
```

For TiraLib-based evaluation pipelines that consume benchmarks from flat zip archives:

```sh
python tools/make_benchmark_archive.py --out polybench_tiramisu.zip [--sizes MINI MEDIUM] [--benchmarks gemm 2mm]
```

The archive contains `<function>/<function>_generator.cpp` + `_init.h` entries plus a `MANIFEST.json` (repo commit, date, variant list) for provenance.

## Benchmarks and sizes

|Benchmark|Description|
|--- |--- |
|2mm		|2 Matrix Multiplications (alpha * A * B * C + beta * D)|
|3mm		|3 Matrix Multiplications ((A.B).(C.D))|
|adi		|Alternating Direction Implicit solver|
|atax		|Matrix Transpose and Vector Multiplication|
|bicg		|BiCG Sub Kernel of BiCGStab Linear Solver|
|cholesky	|Cholesky Decomposition|
|correlation|	Correlation Computation|
|covariance	|Covariance Computation|
|deriche	|	Edge detection filter|
|doitgen	|	Multi-resolution analysis kernel (MADNESS)|
|durbin		|Toeplitz system solver|
|fdtd_2d	|	2-D Finite Different Time Domain Kernel|
|floyd_warshall	|Graph shortest path lengths|
|gemm		|Matrix-multiply C=alpha.A.B+beta.C|
|gemver		|Vector Multiplication and Matrix Addition|
|gesummv		|Scalar, Vector and Matrix Multiplication|
|gramschmidt	|Gram-Schmidt decomposition|
|heat3d		|Heat equation over 3D data domain|
|jacobi1d	|1-D Jacobi stencil computation|
|jacobi2d	|2-D Jacobi stencil computation|
|lu		|LU decomposition|
|ludcmp		|LU decomposition followed by Forward Substitution|
|mvt		|Matrix Vector Product and Transpose|
|nussinov	|Dynamic programming algorithm for sequence alignment|
|seidel2d		|2-D Seidel stencil computation|
|symm		|Symmetric matrix-multiply|
|syr2k		|Symmetric rank-2k operations|
|syrk		|Symmetric rank-k operations|
|trisolv	|	Triangular solver|
|trmm		|Triangular matrix-multiply|

Sizes match PolyBench's datasets exactly: `MINI`, `SMALL`, `MEDIUM`, `LARGE`, `XLARGE` (= PolyBench's `EXTRALARGE`). Data types follow PolyBench: `double` for most benchmarks, `float` for deriche, `int` for floyd_warshall and nussinov.

## Validation

Three validators live in `tools/` (they need a PolyBench/C 4.2.1 checkout — `./tools/get_polybench_c.sh` fetches one):

- `validate_init_parity.py` — proves each init sidecar is numerically identical to stock PolyBench's `init_array`/`print_array` (gcc/g++ only; runs in CI). **Status: 90/90 pass (MINI/SMALL/MEDIUM), 60/60 compile (LARGE/XLARGE).**
- `validate_standalone.py` — full standalone path vs. stock PolyBench/C binaries: init → kernel → live-out dump diff. Requires Tiramisu. **Status: 30/30 match (MINI).**
- `validate_against_polybench_c.py` — same check through TiraLib's `PolybenchHarness`. Requires Tiramisu + TiraLib. **Status: 30/30 match (MINI).**

Known numerical caveat (inherent to the benchmark, not this port): PolyBench's gramschmidt input matrix is rank-deficient, so classical Gram-Schmidt hits a numerically zero pivot (k=13/30 at MINI); beyond it, outputs are ~1e-15 cancellation residue whose value is FP-contraction dependent — stock PolyBench/C itself produces a different tail when compiled with `-mfma -ffp-contract=fast`. The validators compare gramschmidt up to the singular pivot and report the masked region.

## Regenerating the generated files

`*_init.h`, `*_wrapper.cpp` and `*_wrapper.h` are generated from the PolyBench/C sources and the Tiramisu generators:

```sh
python tools/generate_polybench_init_sidecars.py --polybench-c <path/to/PolyBenchC-4.2.1>
```

The generator hard-fails on any buffer↔array mapping, dimension, or type inconsistency. Re-run it (and the validators) after modifying a benchmark.

## References

- Pouchet, Louis-Noël. "Polybench: The polyhedral benchmark suite, 2012." URL: [http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/](http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/) (2012).
- Yuki, Tomofumi, and Louis-Noël Pouchet. "Polybench 4.2.1" (2016).
- Baghdadi, Riyadh, et al. "Tiramisu: A polyhedral compiler for expressing fast and portable code." _2019 IEEE/ACM International Symposium on Code Generation and Optimization (CGO)_. IEEE, 2019.

## License

MIT (see `LICENSE`). `utilities/polybench.{c,h}` and the `init_array`/`print_array` bodies embedded in the generated `*_init.h` files come from PolyBench/C 4.2.1 by Louis-Noël Pouchet and Tomofumi Yuki (see `NOTICE`).
