#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_symm.o.h"
#include "polybench-tiramisu.h"
#include "symm.h"
#include <tiramisu/utils.h>


int symm_ref(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C)
{
  int i,j,k;
  double temp2;
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
        temp2 = 0;
        for (k = 0; k < i; k++) {
            C(k, j) += alpha * B(i, j) * A(i, k);
            temp2 += B(k, j) * A(i, k);
        }
        C(i, j) = beta * C(i, j) + alpha * B(i, j) * A(i, i) + alpha * temp2;
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
    Halide::Buffer<double> b_A(M,M), b_B(N,M), b_C(N,M), b_C_ref(N,M);
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
	    	    symm_ref(b_A, b_B, b_C_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_C_ref);
          transpose(b_A);
          transpose(b_B);
          
          duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_B, b_C);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    symm(b_A.raw_buffer(), b_B.raw_buffer(), b_C.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "symm",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("symm", b_C_ref, b_C, 0.001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_C);
        std::cout << "Reference " << std::endl;
        print_buffer(b_C_ref);
    }

    return 0;
}
