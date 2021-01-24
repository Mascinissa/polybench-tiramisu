#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TSTEPS 20 
#define N 40

#elif TIRAMISU_SMALL

#define TSTEPS 40 
#define N 120

#elif TIRAMISU_MEDIUM

#define TSTEPS 100
#define N 400

#elif TIRAMISU_LARGE

#define TSTEPS 500 
#define N 2000


#elif TIRAMISU_XLARGE

#define TSTEPS 1000 
#define N 4000

#endif 

int init_array(Halide::Buffer<double> A)
{

  transpose(A);

  int i, j, k;
  int n=N;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A(i, j) = ((double) i*(j+2) + 2) / n;

  transpose(A);

  return 0;
}