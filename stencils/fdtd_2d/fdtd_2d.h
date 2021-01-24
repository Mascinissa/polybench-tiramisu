#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TMAX 20 
#define NX 20
#define NY 30

#elif TIRAMISU_SMALL

#define TMAX 40 
#define NX 60
#define NY 80

#elif TIRAMISU_MEDIUM

#define TMAX 100
#define NX 200
#define NY 240

#elif TIRAMISU_LARGE

#define TMAX 500 
#define NX 1000
#define NY 1200

#elif TIRAMISU_XLARGE

#define TMAX 1000 
#define NX 2000
#define NY 2600

#endif 

int init_array(Halide::Buffer<double> ex, Halide::Buffer<double> ey, Halide::Buffer<double> hz, Halide::Buffer<double> fict)
{
  int i, j, k;
  int nx = NX;
  int ny = NY;
  int tmax = TMAX;
  
  transpose(ex);
  transpose(ey);
  transpose(hz);

  for (i = 0; i < tmax; i++)
    fict(i) = (double) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    {
      ex(i, j) = ((double) i*(j+1)) / nx;
      ey(i, j) = ((double) i*(j+2)) / ny;
      hz(i, j) = ((double) i*(j+3)) / nx;
    }

  transpose(ex);
  transpose(ey);
  transpose(hz);

  return 0;
}