#include "polybench-tiramisu.h"

// Problem size
#if TIRAMISU_MINI

#define W 64
#define H 64

#elif TIRAMISU_SMALL

#define W 192
#define H 128

#elif TIRAMISU_MEDIUM

#define W 720
#define H 480

#elif TIRAMISU_LARGE

#define W 4096
#define H 2160

#elif TIRAMISU_XLARGE

#define W 7680
#define H 4320

#endif 

#define alpha 0.25

int init_array(Halide::Buffer<double> imgIn)
{
  int i, j;
  transpose(imgIn);
  //input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < W; i++)
     for (j = 0; j < H; j++)
	imgIn(i, j) = (double) ((313*i+991*j)%65536) / 65535.0f;
  transpose(imgIn);
return 0;
}