#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_syr2k.o.h"
#include "polybench-tiramisu.h"
#include "syr2k.h"
#include <tiramisu/utils.h>


int syr2k_ref(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C)
{
  int i,j,k;
  for (i = 0; i < N; i++){
    for (j = 0; j <= i; j++)
      C(i, j) *= beta;
    for (k = 0; k < M; k++)
      for (j = 0; j <=i ; j++)
        C(i, j) += A(j, k) * alpha*B(i, k) + B(j, k) * alpha*A(i, k);
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
    Halide::Buffer<double> b_A(M,N), b_B(M,N), b_C(N,N), b_C_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_B, b_C_ref);

          transpose(b_C_ref);
          transpose(b_A);
          transpose(b_B);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    syr2k_ref(b_A, b_B, b_C_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_C_ref);
          transpose(b_A);
          transpose(b_B);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_B, b_C);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    syr2k(b_A.raw_buffer(), b_B.raw_buffer(), b_C.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "syr2k",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("syr2k", b_C_ref, b_C, 0.0001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_C);
        std::cout << "Reference " << std::endl;
        print_buffer(b_C_ref);
    }

    return 0;
}
