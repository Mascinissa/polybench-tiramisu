#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define M 28
#define N 32

#elif TIRAMISU_SMALL

#define M 80
#define N 100


#elif TIRAMISU_MEDIUM

#define M 240
#define N 260

#elif TIRAMISU_LARGE

#define M 1200
#define N 1400

#elif TIRAMISU_XLARGE

#define M 2600
#define N 3000

#endif 

int init_array(Halide::Buffer<double> data)
{
  int i, j;

  transpose(data);

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data(i, j) = ((double) i*j) / M;
    
  transpose(data);

return 0;
}