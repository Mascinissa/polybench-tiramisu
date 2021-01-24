#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TSTEPS 20 
#define N 10

#elif TIRAMISU_SMALL

#define TSTEPS 40 
#define N 20

#elif TIRAMISU_MEDIUM

#define TSTEPS 100
#define N 40

#elif TIRAMISU_LARGE

#define TSTEPS 500 
#define N 120


#elif TIRAMISU_XLARGE

#define TSTEPS 1000 
#define N 200

#endif 

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B)
{

  transpose(A);
  transpose(B);
  int i, j, k;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        A(i, j, k) = B(i, j, k) = (double) (i + j + (N-k))* 10 / (N);

  transpose(A);
  transpose(B);
  return 0;
}