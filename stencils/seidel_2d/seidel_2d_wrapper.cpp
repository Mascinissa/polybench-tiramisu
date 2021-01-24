#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_seidel_2d.o.h"
#include "polybench-tiramisu.h"
#include "seidel_2d.h"
#include <tiramisu/utils.h>


int seidel_2d_ref(Halide::Buffer<double> A)
{
int i,j,k,t;
  for (t = 0; t <= TSTEPS - 1; t++)
    for (i = 1; i<= N - 2; i++)
      for (j = 1; j <= N - 2; j++)
	        A(i, j) = (A(i-1, j-1) + A(i-1, j) + A(i-1, j+1)
		                + A(i, j-1) + A(i, j) + A(i, j+1)
		                + A(i+1, j-1) + A(i+1, j) + A(i+1, j+1))/(9.0);  
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
    Halide::Buffer<double> b_A(N,N), b_A_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	        init_array(b_A_ref);
          transpose(b_A_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            seidel_2d_ref(b_A_ref);

          auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        seidel_2d(b_A.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "seidel_2d",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("seidel_2d ", b_A_ref, b_A, 0.00001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
