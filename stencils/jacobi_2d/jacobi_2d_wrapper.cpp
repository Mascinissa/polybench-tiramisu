#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_jacobi_2d.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int jacobi_2d_ref(Halide::Buffer<double> A, Halide::Buffer<double> B)
{
int i,j,k,t;
for (t = 0; t < TSTEPS; t++)
  {
    for (i = 1; i < N - 1; i++)
	    for (j = 1; j < N - 1; j++)
	      B(i, j) = (0.2) * (A(i, j) + A(i, j-1) + A(i, 1+j) + A(1+i, j) + A(i-1, j));
    for (i = 1; i < N - 1; i++)
	    for (j = 1; j < N - 1; j++)
	      A(i, j) = (0.2) * (B(i, j) + B(i, j-1) + B(i, 1+j) + B(1+i, j) + B(i-1, j));
  }   
  return 0;
}

int init_array(Halide::Buffer<double> A, Halide::Buffer<double> B)
{

  transpose(A);
  transpose(B);
  int i, j, k;
  int n=N;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	      A(i, j) = ((double) i*(j+2) + 2) / n;
	      B(i, j) = ((double) i*(j+3) + 3) / n;
      }

  transpose(A);
  transpose(B);
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
    Halide::Buffer<double> b_A(N,N), b_B(N,N), b_A_ref(N,N), b_B_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	        init_array(b_A_ref, b_B_ref);
          transpose(b_A_ref);
          transpose(b_B_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            jacobi_2d_ref(b_A_ref, b_B_ref);

          auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A_ref);
          transpose(b_B_ref);
          
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
	        jacobi_2d(b_A.raw_buffer(), b_B.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "jacobi_2d",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("jacobi_2d A", b_A_ref, b_A, 0.00001);
        compare_buffers_approximately("jacobi_2d B", b_B_ref, b_B, 0.00001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
