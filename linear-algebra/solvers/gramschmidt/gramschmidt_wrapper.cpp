#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_gramschmidt.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int gramschmidt_ref(Halide::Buffer<double> A, Halide::Buffer<double> Q, Halide::Buffer<double> R)
{
  int i,j,k;
  float nrm;
  for (k = 0; k < N; k++)
  {
    nrm = 0.0;
    for (i = 0; i < M; i++)
      nrm += A(i, k) * A(i, k);
    R(k, k) = sqrt(nrm);
    for (i = 0; i < M; i++)
      Q(i, k) = A(i, k) / R(k, k);
    for (j = k + 1; j < N; j++)
    {
      R(k, j) = 0.0;
      for (i = 0; i < M; i++)
        R(k, j) += Q(i, k) * A(i, j);
      for (i = 0; i < M; i++)
        A(i, j) = A(i, j) - Q(i, k) * R(k, j);
    }
  }

  return 0;
}


int init_array(Halide::Buffer<double> A, Halide::Buffer<double> Q, Halide::Buffer<double> R)
{
    int i, j;
  transpose(A);
  transpose(Q);
  transpose(R);
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      A(i, j) = (((double) ((i*j) % M) / M )*100) + 10;
      Q(i, j) = 0.0;
    }
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      R(i, j) = 0.0;

  transpose(A);
  transpose(Q);
  transpose(R);
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
    Halide::Buffer<double> b_A(N,M), b_Q(N,M), b_R(N,N), b_A_ref(N,M), b_Q_ref(N,M), b_R_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < 1; ++i)
        {
	      
          init_array(b_A_ref,b_Q_ref,b_R_ref);
          transpose(b_A_ref);
          transpose(b_Q_ref);
          transpose(b_R_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    gramschmidt_ref(b_A_ref,b_Q_ref,b_R_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A_ref);
          transpose(b_Q_ref);
          transpose(b_R_ref);

          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_Q, b_R);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    gramschmidt(b_A.raw_buffer(), b_Q.raw_buffer(), b_R.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "gramschmidt",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("gramschmidt A", b_A_ref, b_A, 0.001);
        compare_buffers_approximately("gramschmidt R", b_R_ref, b_R, 0.001);
        compare_buffers_approximately("gramschmidt Q", b_Q_ref, b_Q, 0.001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
