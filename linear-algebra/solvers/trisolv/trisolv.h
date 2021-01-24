#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define N 40

#elif TIRAMISU_SMALL

#define N 120

#elif TIRAMISU_MEDIUM

#define N 400

#elif TIRAMISU_LARGE

#define N 2000

#elif TIRAMISU_XLARGE

#define N 4000

#endif 

//initializes a positive semi-definite matrix
int init_array(Halide::Buffer<double> L, Halide::Buffer<double> b, Halide::Buffer<double> x)
{
  int i, j;
  int n = N;

  for (i = 0; i < n; i++)
    {
      x(i) = - 999;
      b(i) =  i ;
      for (j = 0; j <= i; j++)
  	    L(i, j) = (double) (i+n-j+1)*2/n;
    }

return 0;
}