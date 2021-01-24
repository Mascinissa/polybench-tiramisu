#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TSTEPS 20 
#define N 30

#elif TIRAMISU_SMALL

#define TSTEPS 40 
#define N 90

#elif TIRAMISU_MEDIUM

#define TSTEPS 100
#define N 250

#elif TIRAMISU_LARGE

#define TSTEPS 500 
#define N 1300


#elif TIRAMISU_XLARGE

#define TSTEPS 1000 
#define N 2800

#endif 

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B)
{

  transpose(A);
  transpose(B);
  int i, j, k;
  int n=N;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	      A(i, j) = ((double) i*(j+2) + 2) / n;
	      B(i, j) = ((double) i*(j+3) + 3) / n;
      }

  transpose(A);
  transpose(B);
  return 0;
}