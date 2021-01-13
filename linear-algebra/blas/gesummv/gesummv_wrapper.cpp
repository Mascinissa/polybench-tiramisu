#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_gesummv.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int gesummv_ref(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> x, Halide::Buffer<double> y)
{
  int i,j,k;
  Halide::Buffer<double> tmp(N);

    for (i = 0; i < N; i++)
    {
      tmp(i) = 0.0;
      y(i) = 0.0;
      for (j = 0; j < N; j++)
	{
	  tmp(i) = A(i, j) * x(j) + tmp(i);
	  y(i) = B(i, j) * x(j) + y(i);
	}
      y(i) = alpha * tmp(i) + beta * y(i);
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
    Halide::Buffer<double> b_A(N,N), b_B(N,N), b_x(N), b_y(N), b_y_ref(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_buffer(b_y_ref, (double) 5);
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B, (double) 3);
	      init_buffer(b_x, (double) 4);

          transpose(b_A);
          transpose(b_B);
          auto start = std::chrono::high_resolution_clock::now();

	      if (run_ref)
	        gesummv_ref(b_A, b_B, b_x, b_y_ref);

	      auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A);
          transpose(b_B);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_buffer(b_y, (double) 5);
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B, (double) 3);
	      init_buffer(b_x, (double) 4);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        gesummv(b_A.raw_buffer(), b_B.raw_buffer(), b_x.raw_buffer(), b_y.raw_buffer());

	        auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "gesummv",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("gesummv", b_y_ref, b_y);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_y);
        std::cout << "Reference " << std::endl;
        print_buffer(b_y_ref);
    }

    return 0;
}
