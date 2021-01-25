#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_nussinov.o.h"
#include "polybench-tiramisu.h"
#include "nussinov.h"
#include <tiramisu/utils.h>


int nussinov_ref(Halide::Buffer<double> table, Halide::Buffer<double> seq)
{
  int i,j,k;
 for (i = N-1; i >= 0; i--) {
  for (j=i+1; j<N; j++) {


   if (j-1>=0) 
      table(i, j) = std::max(table(i, j), table(i, j-1));
   if (i+1<N) 
      table(i, j) = std::max(table(i, j), table(i+1, j));

   if (j-1>=0 && i+1<N) {
     /* don't allow adjacent elements to bond */
     if (i<j-1) 
        table(i, j) = std::max(table(i, j), table(i+1, j-1)+(double)((seq(i)+seq(j))==3.0));
     else 
        table(i, j) = std::max(table(i, j), table(i+1, j-1));
   }

   for (k=i+1; k<j; k++) {
      table(i, j) = std::max(table(i, j), table(i, k) + table(k+1, j));
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
    Halide::Buffer<double> b_table(N,N), b_table_ref(N,N), b_seq(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_table_ref, b_seq);

          transpose(b_table_ref);
          auto start = std::chrono::high_resolution_clock::now();

	        if (run_ref)
	    	    nussinov_ref(b_table_ref, b_seq);

	        auto end = std::chrono::high_resolution_clock::now();
          transpose(b_table_ref);

          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_table, b_seq);

          auto start = std::chrono::high_resolution_clock::now();
	        if (run_tiramisu)
	    	    nussinov(b_table.raw_buffer(), b_seq.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "nussinov",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("nussinov", b_table_ref, b_table);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_table);
        std::cout << "Reference " << std::endl;
        print_buffer(b_table_ref);
    }

    return 0;
}
