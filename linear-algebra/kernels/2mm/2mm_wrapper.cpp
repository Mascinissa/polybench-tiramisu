#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_2mm.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int b_2mm_ref(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C, Halide::Buffer<double> D)
{
   int i,j,k;
   Halide::Buffer<double> tmp(P,R);

   for (i = 0; i < P; i++)
    for (j = 0; j < R; j++)
      {
	    tmp(i, j) = 0.0;
	    for (k = 0; k < Q; ++k)
	      tmp(i, j) += alpha * A(i, k) * B(k, j);
      }
  for (i = 0; i < P; i++)
    for (j = 0; j < S; j++)
      {
	    D(i, j) *= beta;
	    for (k = 0; k < R; ++k)
	    D(i, j) += tmp(i, k) * C(k, j);
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
    Halide::Buffer<double> b_A(Q,P), b_B(R,Q), b_C(S, R), b_D(S, P), b_D_ref(S, P);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_buffer(b_D_ref, (double) 1);
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B, (double) 3);
	      init_buffer(b_C, (double) 4);

          transpose(b_D_ref);
          transpose(b_A);
          transpose(b_B);
          transpose(b_C);
          auto start = std::chrono::high_resolution_clock::now();

	      if (run_ref)
	        b_2mm_ref(b_A, b_B, b_C, b_D_ref);

	      auto end = std::chrono::high_resolution_clock::now();
          transpose(b_D_ref);
          transpose(b_A);
          transpose(b_B);
          transpose(b_C);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_buffer(b_D, (double) 1);
	      init_buffer(b_A, (double) 2);
	      init_buffer(b_B, (double) 3);
	      init_buffer(b_C, (double) 4);

   
          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        b_2mm(b_A.raw_buffer(), b_B.raw_buffer(), b_C.raw_buffer(), b_D.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "2mm",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("2mm", b_D_ref, b_D);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_D);
        std::cout << "Reference " << std::endl;
        print_buffer(b_D_ref);
    }

    return 0;
}
