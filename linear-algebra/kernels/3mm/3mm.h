#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define P 16
#define R 18
#define Q 20
#define T 22
#define S 24

#elif TIRAMISU_SMALL

#define P 40
#define R 50
#define Q 60
#define T 70
#define S 80

#elif TIRAMISU_MEDIUM

#define P 180
#define R 190
#define Q 200
#define T 210
#define S 220

#elif TIRAMISU_LARGE

#define P 800
#define R 900
#define Q 1000
#define T 1100
#define S 1200

#elif TIRAMISU_XLARGE

#define P 1600
#define R 1800
#define Q 2000
#define T 2200
#define S 2400

#endif 

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C, Halide::Buffer<double> D)
{
  int i, j;
  int ni=P;
  int nk=Q;
  int nj=R;
  int nm=S;
  int nl=T;

  transpose(A);
  transpose(B);
  transpose(C);
  transpose(D);

  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A(i, j) = (double) (i*j % ni) / (5*ni);
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B(i, j) = (double) (i*(j+1) % nj) / (5*nj);
  for (i = 0; i < nj; i++)
    for (j = 0; j < nm; j++)
      C(i, j) = (double) (i*(j+3) % nl) / (5*nl);
  for (i = 0; i < nm; i++)
    for (j = 0; j < nl; j++)
      D(i, j) = (double) (i*(j+2) % nk) / (5*nk);
    
  transpose(A);
  transpose(B);
  transpose(C);
  transpose(D);
return 0;
}