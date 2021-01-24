#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define M 38
#define N 42

#elif TIRAMISU_SMALL

#define M 116
#define N 124

#elif TIRAMISU_MEDIUM

#define M 390
#define N 410

#elif TIRAMISU_LARGE

#define M 1900
#define N 2100

#elif TIRAMISU_XLARGE

#define M 1800
#define N 2200

#endif 

#define alpha 1.5
#define beta 1.2

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> p, Halide::Buffer<double> r)
{
  int i, j;
  int m=M;
  int n=N;
  transpose(A);
 
  for (i = 0; i < m; i++)
    p(i) = (double)(i % m) / m;
  for (i = 0; i < n; i++) {
    r(i) = (double)(i % n) / n;
    for (j = 0; j < m; j++)
      A(i, j) = (double) (i*(j+1) % n)/n;
  }

  transpose(A);
return 0;
}