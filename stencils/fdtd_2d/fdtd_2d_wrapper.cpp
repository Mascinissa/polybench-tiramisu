#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_fdtd_2d.o.h"
#include "polybench-tiramisu.h"
#include "fdtd_2d.h"
#include <tiramisu/utils.h>


int fdtd_2d_ref(Halide::Buffer<double> ex, Halide::Buffer<double> ey, Halide::Buffer<double> hz, Halide::Buffer<double> fict)
{
  int i,j,t;
  for (t = 0; t < TMAX; t++) {
    for (j = 0; j < NY; j++)
      ey(0, j) = fict(t);
    for (i = 1; i < NX; i++)
      for (j = 0; j < NY; j++)
        ey(i, j) = ey(i, j) - 0.5*(hz(i, j) - hz(i - 1, j));
    for (i = 0; i < NX; i++)
      for (j = 1; j < NY; j++)
        ex(i, j) = ex(i, j) - 0.5*(hz(i, j) - hz(i, j - 1));
    for (i = 0; i < NX - 1; i++)
      for (j = 0; j < NY - 1; j++)
        hz(i, j) = hz(i, j) - 0.7*(ex(i, j + 1) - ex(i, j) + ey(i + 1, j) - ey(i, j));
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
    Halide::Buffer<double> b_fict(TMAX), b_ey(NY,NX), b_ex(NY,NX), b_hz(NY,NX), b_ey_ref(NY,NX), b_ex_ref(NY,NX), b_hz_ref(NY,NX);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	        init_array(b_ex_ref, b_ey_ref, b_hz_ref, b_fict);
          transpose(b_ex_ref);
          transpose(b_ey_ref);
          transpose(b_hz_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            fdtd_2d_ref(b_ex_ref, b_ey_ref, b_hz_ref, b_fict);

          auto end = std::chrono::high_resolution_clock::now();
          transpose(b_ex_ref);
          transpose(b_ey_ref);
          transpose(b_hz_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    //TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	        init_array(b_ex, b_ey, b_hz, b_fict);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        fdtd_2d(b_ex.raw_buffer(), b_ey.raw_buffer(), b_hz.raw_buffer(), b_fict.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "fdtd_2d",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("fdtd_2d ex ", b_ex_ref, b_ex, 0.00001);
        compare_buffers_approximately("fdtd_2d ey ", b_ey_ref, b_ey, 0.00001);
        compare_buffers_approximately("fdtd_2d hz ", b_hz_ref, b_hz, 0.00001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_ex);
        std::cout << "Reference " << std::endl;
        print_buffer(b_ex_ref);
    }

    return 0;
}
