#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TSTEPS 20 
#define N 30

#elif TIRAMISU_SMALL

#define TSTEPS 40 
#define N 120

#elif TIRAMISU_MEDIUM

#define TSTEPS 100
#define N 400

#elif TIRAMISU_LARGE

#define TSTEPS 500 
#define N 200


#elif TIRAMISU_XLARGE

#define TSTEPS 1000 
#define N 4000

#endif 

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B)
{

  int i, j, k;
  int n = N;

  for (i = 0; i < n; i++)
    {
	    A(i)  = ((double) i+ 2) / n;
	    B(i)  = ((double) i+ 3) / n;
    }

  return 0;
}