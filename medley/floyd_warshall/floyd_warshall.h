#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define N 60

#elif TIRAMISU_SMALL

#define N 180

#elif TIRAMISU_MEDIUM

#define N 500

#elif TIRAMISU_LARGE

#define N 2800


#elif TIRAMISU_XLARGE

#define N 5600

#endif 

int init_array (Halide::Buffer<double> path)
{
  int i, j;
  transpose(path);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      path(i, j) = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path(i, j) = 999;
    }
  transpose(path);
  return 0;
}