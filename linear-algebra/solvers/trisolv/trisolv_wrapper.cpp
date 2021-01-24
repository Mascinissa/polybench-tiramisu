#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_trisolv.o.h"
#include "polybench-tiramisu.h"
#include "trisolv.h"
#include <tiramisu/utils.h>


int trisolv_ref(Halide::Buffer<double> L, Halide::Buffer<double> b, Halide::Buffer<double> x)
{
  int i,j,k;
  for (i = 0; i < N; i++)
    {
      x(i) = b(i);
      for (j = 0; j <i; j++)
        x(i) -= L(i, j) * x(j);
      x(i) = x(i) / L(i, i);
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
    Halide::Buffer<double> b_L(N,N), b_x_ref(N), b_x(N), b_b(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_L, b_b, b_x_ref);

          transpose(b_L);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    trisolv_ref(b_L, b_b, b_x_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_L);

          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_L, b_b, b_x);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    trisolv(b_L.raw_buffer(), b_b.raw_buffer(), b_x.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "trisolv",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("trisolv", b_x_ref, b_x);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_x);
        std::cout << "Reference " << std::endl;
        print_buffer(b_x_ref);
    }

    return 0;
}
