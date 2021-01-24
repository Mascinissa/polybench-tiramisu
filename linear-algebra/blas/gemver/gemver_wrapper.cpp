#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_gemver.o.h"
#include "polybench-tiramisu.h"
#include "gemver.h"
#include <tiramisu/utils.h>


int gemver_ref(Halide::Buffer<double> A, Halide::Buffer<double> u1, Halide::Buffer<double> u2,
Halide::Buffer<double> v1, Halide::Buffer<double> v2, Halide::Buffer<double> y, Halide::Buffer<double> z,
Halide::Buffer<double> A_hat, Halide::Buffer<double> x, Halide::Buffer<double> w)
{
   int i,j;

   for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      A_hat(i, j) = A(i, j) + u1(i) * v1(j) + u2(i) * v2(j);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      x(i) = x(i) + beta * A_hat(j, i) * y(j);

  for (i = 0; i < N; i++)
    x(i) = x(i) + z(i);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      w(i) = w(i) +  alpha * A_hat(i, j) * x(j);

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
    Halide::Buffer<double> b_A(N, N), b_u1(N), b_u2(N), b_v1(N), b_v2(N), b_z(N), b_y(N), 
                           b_A_hat(N, N), b_A_hat_ref(N, N), b_x(N), b_x_ref(N), b_w(N), b_w_ref(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_u1, b_u2, b_v1, b_v2, b_y, b_z, b_x_ref, b_w_ref);

          transpose(b_A);
          transpose(b_A_hat_ref);
          auto start = std::chrono::high_resolution_clock::now();

	      if (run_ref)
	        gemver_ref(b_A, b_u1, b_u2, b_v1, b_v2, b_y, b_z, b_A_hat_ref, b_x_ref, b_w_ref);

	      auto end = std::chrono::high_resolution_clock::now();
          transpose(b_A);
          transpose(b_A_hat_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_A, b_u1, b_u2, b_v1, b_v2, b_y, b_z, b_x, b_w);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    gemver(b_A.raw_buffer(), b_u1.raw_buffer(), b_u2.raw_buffer(), b_v1.raw_buffer(), 
                b_v2.raw_buffer(), b_y.raw_buffer(), b_z.raw_buffer(), b_A_hat.raw_buffer(), 
                b_x.raw_buffer(), b_w.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }


    print_time("performance_cpu.csv", "gemver",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu){
        compare_buffers_approximately("gemver A_hat", b_A_hat_ref, b_A_hat, 0.001);
        compare_buffers_approximately("gemver x", b_x_ref, b_x, 0.001);
        compare_buffers_approximately("gemver w", b_w_ref, b_w, 0.001);
    }
        

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A_hat);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_hat_ref);
    }

    return 0;
}
