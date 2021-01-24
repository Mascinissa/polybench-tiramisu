#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_heat_3d.o.h"
#include "polybench-tiramisu.h"
#include "heat_3d.h"
#include <tiramisu/utils.h>


int heat_3d_ref(Halide::Buffer<double> A, Halide::Buffer<double> B)
{
  int i,j,k,t;
  for (t = 1; t <= TSTEPS; t++) {
      for (i = 1; i < N-1; i++) {
          for (j = 1; j < N-1; j++) {
              for (k = 1; k < N-1; k++) {
                  B(i, j, k) =   (0.125) * (A(i+1, j, k) - (2.0) * A(i, j, k) + A(i-1, j, k))
                                + (0.125) * (A(i, j+1, k) - (2.0) * A(i, j, k) + A(i, j-1, k))
                                + (0.125) * (A(i, j, k+1) - (2.0) * A(i, j, k) + A(i, j, k-1))
                                + A(i, j, k);
              }
          }
      }
      for (i = 1; i < N-1; i++) {
          for (j = 1; j < N-1; j++) {
              for (k = 1; k < N-1; k++) {
                  A(i, j, k) =   (0.125) * (B(i+1, j, k) - (2.0) * B(i, j, k) + B(i-1, j, k))
                              + (0.125) * (B(i, j+1, k) - (2.0) * B(i, j, k) + B(i, j-1, k))
                              + (0.125) * (B(i, j, k+1) - (2.0) * B(i, j, k) + B(i, j, k-1))
                              + B(i, j, k);
              }
          }
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
    Halide::Buffer<double> b_A(N,N,N), b_B(N,N,N), b_A_ref(N,N,N), b_B_ref(N,N,N);
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
            heat_3d_ref(b_A_ref, b_B_ref);

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
	        heat_3d(b_A.raw_buffer(), b_B.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "heat_3d",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("heat_3d A", b_A_ref, b_A);
        compare_buffers("heat_3d B", b_B_ref, b_B);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
