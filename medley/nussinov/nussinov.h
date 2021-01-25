#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define N 60

#elif TIRAMISU_SMALL

#define N 180

#elif TIRAMISU_MEDIUM

#define N 500

#elif TIRAMISU_LARGE

#define N 2500

#elif TIRAMISU_XLARGE

#define N 5500

#endif 


int init_array(Halide::Buffer<double> table, Halide::Buffer<double> seq)
{
  int i, j;
  int n = N;
  transpose(table);
    //base is AGCT/0..3
  for (i=0; i <n; i++) {
     seq(i) = (char)((i+1)%4);
  }

  for (i=0; i <n; i++)
     for (j=0; j <n; j++)
       table(i, j) = 0;
  transpose(table);

return 0;
}