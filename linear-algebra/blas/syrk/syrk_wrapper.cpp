#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_syrk.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int syrk_ref(Halide::Buffer<double> A, Halide::Buffer<double> C)
{
  int i,j,k;
  for (i = 0; i < N; i++) {
    for (j = 0; j <= i; j++)
      C(i, j) *= beta;
    for (k = 0; k < M; k++) {
      for (j = 0; j <= i; j++)
        C(i, j) += alpha * A(i, k) * A(j, k);
    }
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
    Halide::Buffer<double> b_A(M,N), b_C(N,N), b_C_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_buffer(b_A, (double) 2);
          for (int i = 0; i < M; ++i)  
            for (int j = 0; j < N; ++j)
              b_C_ref(j,i)=pow(-1,(j%2))*i+pow(-1,((j+i)%2))*j;

          transpose(b_C_ref);
          transpose(b_A);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    syrk_ref(b_A, b_C_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_C_ref);
          transpose(b_A);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_buffer(b_A, (double) 2);
          for (int i = 0; i < M; ++i)  
            for (int j = 0; j < N; ++j)
              b_C(j,i)=pow(-1,(j%2))*i+pow(-1,((j+i)%2))*j;

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    syrk(b_A.raw_buffer(), b_C.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "syrk",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("syrk", b_C_ref, b_C);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_C);
        std::cout << "Reference " << std::endl;
        print_buffer(b_C_ref);
    }

    return 0;
}
