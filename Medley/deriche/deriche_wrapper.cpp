#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_deriche.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>

int deriche_ref(Halide::Buffer<double> imgIn, Halide::Buffer<double> imgOut)
{
  int i,j;
  Halide::Buffer<double> xm1(W), tm1(W), ym1(W), ym2(W);
  Halide::Buffer<double> xp1(W), xp2(W);
  Halide::Buffer<double> tp1(W), tp2(W);
  Halide::Buffer<double> yp1(W), yp2(W);
  Halide::Buffer<double> y1(W,H), y2(W,H);
  
  double k;
  double a1, a2, a3, a4, a5, a6, a7, a8;
  double b1, b2, c1, c2;
  k = (1.0-exp(-alpha))*(1.0-exp(-alpha))/(1.0+2.0*alpha*exp(-alpha)-exp(2.0*alpha));
  a1 = a5 = k;
  a2 = a6 = k*exp(-alpha)*(alpha-1.0);
  a3 = a7 = k*exp(-alpha)*(alpha+1.0);
  a4 = a8 = -k*exp((-2.0)*alpha);
  b1 =  pow(2.0,-alpha);
  b2 = -exp((-2.0)*alpha);
  c1 = c2 = 1;

  for (i=0; i<W; i++) {
      ym1(i) = 0.0;
      ym2(i) = 0.0;
      xm1(i) = 0.0;
      for (j=0; j<H; j++) {
          y1(i, j) = a1*imgIn(i, j) + a2*xm1(i) + b1*ym1(i) + b2*ym2(i);
          xm1(i) = imgIn(i, j);
          ym2(i) = ym1(i);
          ym1(i) = y1(i, j);
      }
  }

  for (i=0; i<W; i++) {
      yp1(i) = 0.0;
      yp2(i) = 0.0;
      xp1(i) = 0.0;
      xp2(i) = 0.0;
      for (j=0; j<H; j++) {
          y2(i, H-1-j) = a3*xp1(i) + a4*xp2(i) + b1*yp1(i) + b2*yp2(i);
          xp2(i) = xp1(i);
          xp1(i) = imgIn(i, H-1-j);
          yp2(i) = yp1(i);
          yp1(i) = y2(i, H-1-j);
      }
  }

  for (i=0; i<W; i++)
      for (j=0; j<H; j++) {
          imgOut(i, j) = c1 * (y1(i, j) + y2(i, j));
      }

  for (j=0; j<H; j++) {
      tm1(j) = 0.0;
      ym1(j) = 0.0;
      ym2(j) = 0.0;
      for (i=0; i<W; i++) {
          y1(i, j) = a5*imgOut(i, j) + a6*tm1(j) + b1*ym1(j) + b2*ym2(j);
          tm1(j) = imgOut(i, j);
          ym2(j) = ym1(j);
          ym1(j) = y1 (i, j);
      }
  }
  

  for (j=0; j<H; j++) {
      tp1(j) = 0.0;
      tp2(j) = 0.0;
      yp1(j) = 0.0;
      yp2(j) = 0.0;
      for (i=0; i<W; i++) {
          y2(W-1-i, j) = a7*tp1(j) + a8*tp2(j) + b1*yp1(j) + b2*yp2(j);
          tp2(j) = tp1(j);
          tp1(j) = imgOut(W-1-i, j);
          yp2(j) = yp1(j);
          yp1(j) = y2(W-1-i, j);
      }
  }

  for (i=0; i<W; i++)
      for (j=0; j<H; j++)
          imgOut(i, j) = c2*(y1(i, j) + y2(i, j));
  return 0;
}

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

int main(int argc, char** argv)
{
    std::vector<std::chrono::duration<double, std::milli>> duration_vector_1, duration_vector_2;
    bool run_ref = false, run_tiramisu = false;

    const char* env_ref = std::getenv("RUN_REF");

    if (env_ref != NULL && env_ref[0] == '1')
        run_ref = true;

    const char* env_tiramisu = std::getenv("RUN_TIRAMISU");

    if (env_tiramisu != NULL && env_tiramisu[0] == '1')
        run_tiramisu = true;

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    Halide::Buffer<double> b_imgIn(H,W), b_imgOut_ref(H,W), b_imgOut(H,W);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_imgIn);

          transpose(b_imgIn);
          transpose(b_imgOut_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    deriche_ref(b_imgIn, b_imgOut_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_imgIn);
          transpose(b_imgOut_ref);

          duration_vector_1.push_back(end - start);
        }
    }
        
    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_imgIn);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    deriche(b_imgIn.raw_buffer(), b_imgOut.raw_buffer());

  	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
          
        }
    }

    print_time("performance_cpu.csv", "deriche",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("deriche", b_imgOut_ref, b_imgOut, 0.0001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_imgOut);
        std::cout << "Reference " << std::endl;
        print_buffer(b_imgOut_ref);
    }

    return 0;
}
