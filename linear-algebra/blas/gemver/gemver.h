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

#define alpha 1.5
#define beta 1.2



int init_array(Halide::Buffer<double> A, Halide::Buffer<double> u1, Halide::Buffer<double> u2,
Halide::Buffer<double> v1, Halide::Buffer<double> v2, Halide::Buffer<double> y, Halide::Buffer<double> z,
Halide::Buffer<double> x, Halide::Buffer<double> w)
{
  int i, j;
  double fn = (double) N;
  int n = N;
  transpose(A);
    for (i = 0; i < n; i++)
    {   
      u1(i) = i;
      u2(i) = ((i+1)/fn)/2.0;
      v1(i) = ((i+1)/fn)/4.0;
      v2(i) = ((i+1)/fn)/6.0;
      y(i) = ((i+1)/fn)/8.0;
      z(i) = ((i+1)/fn)/9.0;
      x(i) = 0.0;
      w(i) = 0.0;
      for (j = 0; j < n; j++)
        A(i, j) = (double) (i*j % n) / n;
    }   
  transpose(A);
return 0;
}