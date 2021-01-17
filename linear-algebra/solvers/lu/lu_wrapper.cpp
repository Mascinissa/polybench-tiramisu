#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_lu.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int lu_ref(Halide::Buffer<double> A)
{
  int i,j,k;
  for (i = 0; i < N; i++) {
    for (j = 0; j <i; j++) {
       for (k = 0; k < j; k++) {
          A(i, j) -= A(i, k) * A(k, j);
       }
        A(i, j) /= A(j, j);
    }
   for (j = i; j < N; j++) {
       for (k = 0; k < i; k++) {
          A(i, j) -= A(i, k) * A(k, j);
       }
    }
  }
  return 0;
}

//initializes a positive semi-definite matrix
int init_array(Halide::Buffer<double> A)
{
  int i, j;
  int n = N;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A(i, j) = (double)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A(i, j) = 0;
      }
      A(i, i) = 1;
    }

  /* Make the matrix positive semi-definite. */
  int r,s,t;
  Halide::Buffer<double> B(N,N);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      B(r, s) = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	B(r, s) += A(r, t) * A(s, t);
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A(r, s) = B(r, s);

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
	    	    lu_ref(b_A_ref);

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
	    	    lu(b_A.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "lu",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("lu", b_A_ref, b_A);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
