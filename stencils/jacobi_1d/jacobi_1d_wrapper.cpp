#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_jacobi_1d.o.h"
#include "polybench-tiramisu.h"
#include "jacobi_1d.h"
#include <tiramisu/utils.h>


int jacobi_1d_ref(Halide::Buffer<double> A, Halide::Buffer<double> B)
{
  int i,j,k,t;
  for (t = 0; t < TSTEPS; t++)
    {
      for (i = 1; i < N - 1; i++)
	      B(i)  = 0.33333 * (A(i-1)  + A(i)  + A(i + 1) );
      for (i = 1; i < N - 1; i++)
	      A(i)  = 0.33333 * (B(i-1)  + B(i)  + B(i + 1) );
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
    Halide::Buffer<double> b_A(N), b_B(N), b_A_ref(N), b_B_ref(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	        init_array(b_A_ref, b_B_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            jacobi_1d_ref(b_A_ref, b_B_ref);

          auto end = std::chrono::high_resolution_clock::now();

          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_B);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        jacobi_1d(b_A.raw_buffer(), b_B.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "jacobi_1d",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("jacobi_1d A", b_A_ref, b_A, 0.00001);
        compare_buffers_approximately("jacobi_1d B", b_B_ref, b_B, 0.00001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
