#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_covariance.o.h"
#include "polybench-tiramisu.h"
#include "covariance.h"
#include <tiramisu/utils.h>

//TODO
int covariance_ref(Halide::Buffer<double> data, Halide::Buffer<double> cov)
{
    Halide::Buffer<double> mean(M);
    int i,j,k;

    for (j = 0; j < M; j++)
    {
      mean(j) = 0.0;
      for (i = 0; i < N; i++)
        mean(j) += data(i, j);
      mean(j) /= N;
    }

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data(i, j) -= mean(j);

  for (i = 0; i < M; i++)
    for (j = i; j < M; j++)
      {
        cov(i, j) = 0.0;
        for (k = 0; k < N; k++)
	  cov(i, j) += data(k, i) * data(k, j);
        cov(i, j) /= (N - 1.0);
        cov(j, i) = cov(i, j);
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
    Halide::Buffer<double> b_data(M,N), b_cov(M,M), b_cov_ref(M, M);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	        init_array(b_data);
          
          transpose(b_data);
          transpose(b_cov_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    covariance_ref(b_data, b_cov_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_data);
          transpose(b_cov_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	        init_array(b_data);

            auto start = std::chrono::high_resolution_clock::now();

	        if (run_tiramisu)
	    	    covariance(b_data.raw_buffer(), b_cov.raw_buffer());

	        auto end = std::chrono::high_resolution_clock::now();
            duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "covariance",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("covariance", b_cov_ref, b_cov, 0.0001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_cov);
        std::cout << "Reference " << std::endl;
        print_buffer(b_cov_ref);
    }

    return 0;
}
