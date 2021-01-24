#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_3mm.o.h"
#include "polybench-tiramisu.h"
#include "3mm.h"
#include <tiramisu/utils.h>


int b_3mm_ref(Halide::Buffer<double> A, Halide::Buffer<double> B, Halide::Buffer<double> C, Halide::Buffer<double> D, Halide::Buffer<double> G)
{
   int i,j,k;
   Halide::Buffer<double> E(P,R), F(R,T);

  /* E := A*B */
  for (i = 0; i < P; i++)
    for (j = 0; j < R; j++)
      {
	E(i, j) = 0.0;
	for (k = 0; k < Q; ++k)
	  E(i, j) += A(i, k) * B(k, j);
      }
  /* F := C*D */
  for (i = 0; i < R; i++)
    for (j = 0; j < T; j++)
      {
	F(i, j) = 0.0;
	for (k = 0; k < S; ++k)
	  F(i, j) += C(i, k) * D(k, j);
      }
  /* G := E*F */
  for (i = 0; i < P; i++)
    for (j = 0; j < T; j++)
      {
	G(i, j) = 0.0;
	for (k = 0; k < R; ++k)
	  G(i, j) += E(i, k) * F(k, j);
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
    Halide::Buffer<double> b_A(Q,P), b_B(R,Q), b_C(S, R), b_D(T, S), b_E_ref(T, P), b_E(T, P);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_array(b_A, b_B, b_C, b_D);

          transpose(b_E_ref);
          transpose(b_A);
          transpose(b_B);
          transpose(b_C);
          transpose(b_D);
          auto start = std::chrono::high_resolution_clock::now();

	      if (run_ref)
	        b_3mm_ref(b_A, b_B, b_C, b_D, b_E_ref);

	      auto end = std::chrono::high_resolution_clock::now();
          transpose(b_E_ref);
          transpose(b_A);
          transpose(b_B);
          transpose(b_C);
          transpose(b_D);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_array(b_A, b_B, b_C, b_D);
   
          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	          b_3mm(b_A.raw_buffer(), b_B.raw_buffer(), b_C.raw_buffer(), b_D.raw_buffer(), b_E.raw_buffer());

	        auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "3mm",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("3mm", b_E_ref, b_E, 0.0001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_E);
        std::cout << "Reference " << std::endl;
        print_buffer(b_E_ref);
    }

    return 0;
}
