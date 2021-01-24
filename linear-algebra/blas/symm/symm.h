#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define M 5
#define N 6

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

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C)
{
    int i, j;
    int m=M;
    int n=N;
  transpose(A);
  transpose(B);
  transpose(C);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      C(i, j) = (double) ((i+j) % 100) / m;
      B(i, j) = (double) ((n+i-j) % 100) / m;
    }
  for (i = 0; i < m; i++) {
    for (j = 0; j <=i; j++)
      A(i, j) = (double) ((i+j) % 100) / m;
    for (j = i+1; j < m; j++)
      A(i, j) = -999; //regions of arrays that should not be used
  }

  transpose(A);
  transpose(B);
  transpose(C);
return 0;
}