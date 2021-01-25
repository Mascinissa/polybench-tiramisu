#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_durbin.o.h"
#include "polybench-tiramisu.h"
#include "durbin.h"
#include <tiramisu/utils.h>


int durbin_ref(Halide::Buffer<double> r, Halide::Buffer<double> y)
{
  int i,j,k;
  double alpha, beta, sum;
  Halide::Buffer<double> z(N);
  y(0) = -r(0);
  beta = (1.0);
  alpha= -r(0);

 for (k = 1; k < N; k++) {
   beta = (1-alpha*alpha)*beta;
   sum = (0.0);
   for (i=0; i<k; i++) {
      sum += r(k-i-1)*y(i);
   }
   alpha = - (r(k) + sum)/beta;

   for (i=0; i<k; i++) {
      z(i) = y(i) + alpha*y(k-i-1);
   }
   for (i=0; i<k; i++) {
     y(i) = z(i);
   }
   y(k) = alpha;
 }
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
    Halide::Buffer<double> b_y(N), b_y_ref(N), b_r(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_r);

          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    durbin_ref(b_r, b_y_ref);

	        auto end = std::chrono::high_resolution_clock::now();

          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_r);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    durbin(b_r.raw_buffer(), b_y.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }

    }

    print_time("performance_cpu.csv", "durbin",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("durbin", b_y_ref, b_y, 0.0001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_y);
        std::cout << "Reference " << std::endl;
        print_buffer(b_y_ref);
    }

    return 0;
}
