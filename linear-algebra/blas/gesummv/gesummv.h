#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define N 30

#elif TIRAMISU_SMALL

#define N 90

#elif TIRAMISU_MEDIUM

#define N 250

#elif TIRAMISU_LARGE

#define N 1300

#elif TIRAMISU_XLARGE

#define N 2800

#endif 

#define alpha 1.5
#define beta 1.2


int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> x)
{
  int i, j;
  int n = N;
  transpose(A);
  transpose(B);

    for (i = 0; i < n; i++)
    {
      x(i) = (double)( i % n) / n;
      for (j = 0; j < n; j++) {
	      A(i, j) = (double) (i*j % n) / n;
	      B(i, j) = (double) (i*j % n) / n;
      }
    }
  transpose(A);
  transpose(B);
return 0;
}