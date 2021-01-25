#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_correlation.o.h"
#include "polybench-tiramisu.h"
#include "correlation.h"
#include <tiramisu/utils.h>


int correlation_ref(Halide::Buffer<double> data, Halide::Buffer<double> corr)
{
    Halide::Buffer<double> mean(M),stddev(M);
    int i,j,k;
    // double  eps = 0.1;

   for (j = 0; j < M; j++)
    {
      mean(j) = (0.0);
      for (i = 0; i < N; i++)
	    mean(j) += data(i, j);
      mean(j) /= N;
    }


   for (j = 0; j < M; j++)
    {
      stddev(j) = (0.0);
      for (i = 0; i < N; i++)
        stddev(j) += (data(i, j) - mean(j)) * (data(i, j) - mean(j));
      stddev(j) /= N;
      stddev(j) = sqrt(stddev(j));
      /* The following in an inelegant but usual way to handle
         near-zero std. dev. values, which below would cause a zero-
         divide. */
      stddev(j) = stddev(j) <= eps ? (1.0) : stddev(j);
    }

  /* Center and reduce the column vectors. */
  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      {
        data(i, j) -= mean(j);
        data(i, j) /= sqrt(N) * stddev(j);
      }

  /* Calculate the m * m correlation matrix. */
  for (i = 0; i < M-1; i++)
    {
      corr(i, i) = (1.0);
      for (j = i+1; j < M; j++)
        {
          corr(i, j) = (0.0);
          for (k = 0; k < N; k++)
            corr(i, j) += (data(k, i) * data(k, j));
          corr(j, i) = corr(i, j);
        }
    }
  corr(M-1, M-1) = (1.0);

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
    Halide::Buffer<double> b_data(M,N), b_corr(M,M), b_corr_ref(M, M);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
          init_array(b_data);
          
          transpose(b_data);
          transpose(b_corr_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    correlation_ref(b_data, b_corr_ref);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_data);
          transpose(b_corr_ref);
          
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
	    	    correlation(b_data.raw_buffer(), b_corr.raw_buffer());

	        auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "correlation",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("correlation", b_corr_ref, b_corr, 0.001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_corr);
        std::cout << "Reference " << std::endl;
        print_buffer(b_corr_ref);
    }

    return 0;
}
