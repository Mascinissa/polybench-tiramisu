#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_trmm.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int trmm_ref(Halide::Buffer<double> A, Halide::Buffer<double> B)
{
  int i,j,k;
  for (i = 0; i < M; i++)
     for (j = 0; j < N; j++) {
        for (k = i+1; k < M; k++) 
           B(i, j) += A(k, i) * B(k, j);
        B(i, j) = alpha * B(i, j);
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
    Halide::Buffer<double> b_A(M,M), b_B(N,M), b_B_ref(N,M);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B_ref, (double) 3);


          transpose(b_B_ref);
          transpose(b_A);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    trmm_ref(b_A, b_B_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_B_ref);
          transpose(b_A);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B, (double) 3);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    trmm(b_A.raw_buffer(), b_B.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "trmm",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("trmm", b_B_ref, b_B);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_B);
        std::cout << "Reference " << std::endl;
        print_buffer(b_B_ref);
    }

    return 0;
}
