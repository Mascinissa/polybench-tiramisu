# polybench-tiramisu
 PolyBench-Tiramisu is the reimplementation of [PolyBench](http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/) 4.2  in the [Tiramisu](http://tiramisu-compiler.org/) programming language. It is a benchmark suite of 30 numerical computations with static control flow, extracted from operations in various application domains (linear algebra computations, image processing, physics simulation, dynamic programming, statistics, etc.).
 
## Usage 
Set the environment variable `TIRAMISU_ROOT` to the root directory of your Tiramisu installation.

    export TIRAMISU_ROOT=<path/to/Tiramisu>
Run compile_and_run.sh specifying the benchmark name and the problem size.

    ./compile_and_run.sh <benchmark name> <problem size>

Example:

    ./compile_and_run.sh cholesky LARGE

### Available Benchmarks
|Benchmark|Description|
|--- |--- |
|2mm		|2 Matrix Multiplications (E=A.B; F=E.C; G=F+D)|
|3mm		|3 Matrix Multiplications (E=A.B; F=C.D; G=E.F)|
|adi		|Alternating Direction Implicit solver|
|atax		|Matrix Transpose and Vector Multiplication|
|bicg		|BiCG Sub Kernel of BiCGStab Linear Solver|
|cholesky	|Cholesky Decomposition|
|correlation|	Correlation Computation|
|covariance	|Covariance Computation|
|deriche	|	Edge detection filter|
|doitgen	|	Multi-resolution analysis kernel (MADNESS)|
|durbin		|Toeplitz system solver|
|fdtd-2d	|	2-D Finite Different Time Domain Kernel|
|gemm		|Matrix-multiply C=alpha.A.B+beta.C|
|gemver		|Vector Multiplication and Matrix Addition|
|gesummv		|Scalar, Vector and Matrix Multiplication|
|gramschmidt	|Gram-Schmidt decomposition|
|head-3d		|Heat equation over 3D data domain|
|jacobi-1D	|1-D Jacobi stencil computation|
|jacobi-2D	|2-D Jacobi stencil computation|
|lu		|LU decomposition|
|ludcmp		|LU decomposition followed by FS|
|mvt		|Matrix Vector Product and Transpose|
|nussinov	|Dynamic programming algorithm for sequence alignment|
|seidel		|2-D Seidel stencil computation|
|symm		|Symmetric matrix-multiply|
|syr2k		|Symmetric rank-2k operations|
|syrk		|Symmetric rank-k operations|
|trisolv	|	Triangular solver|
|trmm		|Triangular matrix-multiply|

### Available Sizes
| Problem size | Memory usage |
|--|--|
| MINI | < 16KB of memory. |
| SMALL| ≈ 128KB of memory.|
| MEDIUM| ≈ 1MB of memory.|
| LARGE| ≈ 25MB of memory.|
| EXTRALARGE| ≈ 120MB of memory.|


##	 References

 - Pouchet, Louis-Noël. "Polybench: The polyhedral benchmark suite,
   2012." URL: [http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/](http://web.cse.ohio-state.edu/~pouchet.2/software/polybench/)  (2012).
  - Yuki, Tomofumi, and Louis-Noël Pouchet. "Polybench 4.0." (2015).
   
- Baghdadi, Riyadh, et al. "Tiramisu: A polyhedral compiler for expressing fast and portable code." _2019 IEEE/ACM International Symposium on Code Generation and Optimization (CGO)_. IEEE, 2019.
