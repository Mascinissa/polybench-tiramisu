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


int init_array(Halide::Buffer<double> A, Halide::Buffer<double> y1, Halide::Buffer<double> y2, Halide::Buffer<double> x1, Halide::Buffer<double> x2)
{
  int i, j;
  int n = N;
  transpose(A);
  for (i = 0; i < n; i++)
    {
      x1(i) = (double) (i % n) / n;
      x2(i) = (double) ((i + 1) % n) / n;
      y1(i) = (double) ((i + 3) % n) / n;
      y2(i) = (double) ((i + 4) % n) / n;
      for (j = 0; j < n; j++)
	      A(i, j) = (double) (i*j % n) / n;
    }
  transpose(A);
return 0;
}