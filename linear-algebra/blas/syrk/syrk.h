#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define M 20
#define N 30

#elif TIRAMISU_SMALL

#define M 60
#define N 80

#elif TIRAMISU_MEDIUM

#define M 200
#define N 240

#elif TIRAMISU_LARGE

#define M 1000
#define N 1200

#elif TIRAMISU_XLARGE

#define M 2000
#define N 2600

#endif 

#define alpha 1.5
#define beta 1.2

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> C)
{
    int i, j;
    int m=M;
    int n=N;
  transpose(A);
  transpose(C);

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      A(i, j) = (double) (i*j%n) / n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C(i, j) = (double) (i*j%m) / m;

  transpose(A);
  transpose(C);
return 0;
}