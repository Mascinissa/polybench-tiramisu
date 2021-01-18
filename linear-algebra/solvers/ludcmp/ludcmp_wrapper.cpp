#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_ludcmp.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int ludcmp_ref(Halide::Buffer<double> A, Halide::Buffer<double> b, Halide::Buffer<double> y, Halide::Buffer<double> x)
{
  int i,j,k;
  double w;
    for (i = 0; i < N; i++) {
    for (j = 0; j <i; j++) {
       w = A(i, j);
       for (k = 0; k < j; k++) {
          w -= A(i, k) * A(k, j);
       }
        A(i, j) = w / A(j, j);
    }
   for (j = i; j < N; j++) {
       w = A(i, j);
       for (k = 0; k < i; k++) {
          w -= A(i, k) * A(k, j);
       }
       A(i, j) = w;
    }
  }

  for (i = 0; i < N; i++) {
     w = b(i);
     for (j = 0; j < i; j++)
        w -= A(i, j) * y(j);
     y(i) = w;
  }

   for (i = N-1; i >=0; i--) {
     w = y(i);
     for (j = i+1; j < N; j++)
        w -= A(i, j) * x(j);
     x(i) = w / A(i, i);
  }
  return 0;
}

//initializes a positive semi-definite matrix
int init_array(Halide::Buffer<double> A, Halide::Buffer<double> b, Halide::Buffer<double> y, Halide::Buffer<double> x)
{
  int i, j;
  int n = N;
  double fn = (double)n;

  for (i = 0; i < n; i++)
    {
      x(i) = 0;
      y(i) = 0;
      b(i) = (i+1)/fn/2.0 + 4;
    }

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
  /* not necessary for LU, but using same code as cholesky */
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
    Halide::Buffer<double> b_b(N), b_y(N), b_x(N), b_x_ref(N), b_A(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          transpose(b_A);
          init_array(b_A, b_b, b_y, b_x_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    ludcmp_ref(b_A, b_b, b_y, b_x_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A);

          duration_vector_1.push_back(end - start);
        }
    }

   // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_array(b_A, b_b, b_y, b_x);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    ludcmp(b_A.raw_buffer(), b_b.raw_buffer(), b_y.raw_buffer(), b_x.raw_buffer());

	        auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "ludcmp",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("ludcmp", b_x_ref, b_x);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_x);
        std::cout << "Reference " << std::endl;
        print_buffer(b_x_ref);
    }

    return 0;
}
