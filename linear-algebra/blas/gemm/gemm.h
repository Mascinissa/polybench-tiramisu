#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define P 20
#define R 25
#define Q 30

#elif TIRAMISU_SMALL

#define P 60
#define R 70
#define Q 80

#elif TIRAMISU_MEDIUM

#define P 200
#define R 220
#define Q 240

#elif TIRAMISU_LARGE

#define P 1000
#define R 1100
#define Q 1200

#elif TIRAMISU_XLARGE

#define P 2000
#define R 2300
#define Q 2600

#endif 

#define alpha 1.5
#define beta 1.2



int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C)
{
  int i, j;
  int ni = P;
  int nj = R;
  int nk = Q;
  transpose(A);
  transpose(B);
  transpose(C);

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C(i, j) = (double) (i*j % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A(i, j) = (double) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B(i, j) = (double) (i*(j+2) % nj) / nj;

  transpose(A);
  transpose(B);
  transpose(C);
return 0;
}