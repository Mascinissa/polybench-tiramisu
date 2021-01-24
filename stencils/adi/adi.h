#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define TSTEPS 20 
#define N 20

#elif TIRAMISU_SMALL

#define TSTEPS 40 
#define N 60

#elif TIRAMISU_MEDIUM

#define TSTEPS 100
#define N 200

#elif TIRAMISU_LARGE

#define TSTEPS 500 
#define N 1000


#elif TIRAMISU_XLARGE

#define TSTEPS 1000 
#define N 2000

#endif 

int init_array(Halide::Buffer<double> u)
{
  int i, j;
transpose(u);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
    {
	    u(i, j) =  (double)(i + N-j) / N;
    }
transpose(u);
  return 0;
}