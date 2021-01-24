#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define Q 8
#define R 10
#define P 12

#elif TIRAMISU_SMALL

#define Q 20
#define R 25
#define P 30

#elif TIRAMISU_MEDIUM

#define Q 40
#define R 50
#define P 60

#elif TIRAMISU_LARGE

#define Q 140
#define R 150
#define P 160

#elif TIRAMISU_XLARGE

#define Q 220
#define R 250
#define P 270

#endif 


int init_array(Halide::Buffer<double> A, Halide::Buffer<double> x)
{
  int i, j, k;
  int nr = R;
  int np = P;
  int nq = Q;
  transpose(A);
  transpose(x);
   for (i = 0; i < nr; i++)
    for (j = 0; j < nq; j++)
      for (k = 0; k < np; k++)
	      A(i, j, k) = (double) ((i*j + k)%np) / np;
  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++)
      x(i, j) = (double) (i*j % np) / np;
  transpose(A);
  transpose(x);
return 0;
}