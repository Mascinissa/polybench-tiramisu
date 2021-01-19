#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_floyd_warshall.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>


int floyd_warshall_ref(Halide::Buffer<double> path)
{
  int i,j,k;

  for (k = 0; k < N; k++)
    {
      for(i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            path(i, j) = std::min(path(i, j), path(i, k) + path(k, j));
    }
  return 0;
}

int init_array (Halide::Buffer<double> path)
{
  int i, j;
  transpose(path);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      path(i, j) = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path(i, j) = 999;
    }
  transpose(path);
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
    Halide::Buffer<double> b_paths(N,N), b_paths_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	      init_array(b_paths_ref);
          transpose(b_paths_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            floyd_warshall_ref(b_paths_ref);

          auto end = std::chrono::high_resolution_clock::now();
          transpose(b_paths_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_paths);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        floyd_warshall(b_paths.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "floyd_warshall",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("floyd_warshall", b_paths_ref, b_paths);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_paths);
        std::cout << "Reference " << std::endl;
        print_buffer(b_paths_ref);
    }

    return 0;
}
