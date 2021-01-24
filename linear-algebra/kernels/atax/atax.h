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

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> x)
{
  int i, j;
  int m=M;
  int n=N;
  double fn = (double) n;

  transpose(A);

  for (i = 0; i < n; i++)
      x(i) = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A(i, j) = (double) ((i+j) % n) / (5*m);
    
  transpose(A);
return 0;
}