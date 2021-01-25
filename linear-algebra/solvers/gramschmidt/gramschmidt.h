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

#define eps 0.000001

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> Q, Halide::Buffer<double> R)
{
    int i, j;
  transpose(A);
  transpose(Q);
  transpose(R);
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      A(i, j) = (((double) ((i*j) % M) / M )*100) + 10;
      Q(i, j) = 0.0;
    }
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      R(i, j) = 0.0;

  transpose(A);
  transpose(Q);
  transpose(R);
return 0;
}