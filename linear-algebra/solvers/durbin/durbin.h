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


int init_array(Halide::Buffer<double> r)
{
 int i, j;
 int n = N;
 for (i = 0; i < N; i++)
    {
      r(i) = (N+1-i);
    }

return 0;
}