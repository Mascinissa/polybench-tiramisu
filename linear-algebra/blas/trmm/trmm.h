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

#define alpha 1.5


int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B)
{
    int i, j;
  transpose(A);
  transpose(B);

  for (i = 0; i < M; i++) {
    for (j = 0; j < i; j++) {
      A(i, j) = (double)((i+j) % M)/M;
    }
    A(i, i) = 1.0;
    for (j = 0; j < N; j++) {
      B(i, j) = (double)((N+(i-j)) % N)/N;
    }
 }

  transpose(A);
  transpose(B);
return 0;
}