#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_adi.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>

int adi_ref(Halide::Buffer<double> u)
{
  int t, i, j;
  double DX, DY, DT;
  double B1, B2;
  double mul1, mul2;
  double a, b, c, d, e, f;
  Halide::Buffer<double> v(N,N); 
  Halide::Buffer<double> p(N,N); 
  Halide::Buffer<double> q(N,N);

  DX = 1.0/(double)N;
  DY = 1.0/(double)N;
  DT = 1.0/(double)TSTEPS;
  B1 = 2.0;
  B2 = 1.0;
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 /  2.0;
  b = 1.0+mul1;
  c = a;
  d = -mul2 / 2.0;
  e = 1.0+mul2;
  f = d;

for (i=0; i<N; i++) {
  for (j=0; j<N; j++) {
    v(i,j)=1.0;
    q(i,j)=1.0;
    p(i,j)=0.0;
  }
}

 for (t=1; t<=TSTEPS; t++) {
    //Column Sweep
    for (i=1; i<N-1; i++) {
      for (j=1; j<N-1; j++) {
        p(i, j) = -c / (p(i, j-1)*a+b);
        q(i, j) = ((u(j, i-1)*(-d)+u(j, i)*(1.0+2.0*d) - u(j, i+1)*f-q(i, j-1)*a)/(p(i, j-1)*a+b));
      }
      
      for (j=1; j<N-1; j++) {
        v(i, j) = p(i, j) * v(i, j+1) + q(i, j);
      }
    }
    //Row Sweep
    for (i=1; i<N-1; i++) {
      for (j=1; j<N-1; j++) {
        p(i, j) = -f / (p(i, j-1)*d+e);
        q(i, j) = (v(j, i-1)*(-a)+v(j, i)*(1.0+2.0*a) - v(j, i+1)*c-q(i, j-1)*d)/(p(i, j-1)*d+e);
      }
      for (j=1; j<N-1; j++) {
        u(i, j) = p(i, j) * u(i, j+1) + q(i, j);
      }
    }
  }
  return 0;
}

int init_array(Halide::Buffer<double> u)
{
  int i, j;
transpose(u);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
    {
	    u(i, j) =  (double)(i + N-j) / N;
    }
transpose(u);
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
    Halide::Buffer<double> b_u(N,N), b_u_ref(N,N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    //REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {

	      init_array(b_u_ref);
          transpose(b_u_ref);
          auto start = std::chrono::high_resolution_clock::now();

          if (run_ref)
            adi_ref(b_u_ref);

          auto end = std::chrono::high_resolution_clock::now();
          transpose(b_u_ref);
          
          duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i)
        {
	      init_array(b_u);

          auto start = std::chrono::high_resolution_clock::now();
	      if (run_tiramisu)
	        adi(b_u.raw_buffer());

	      auto end = std::chrono::high_resolution_clock::now();
          duration_vector_2.push_back(end - start);
        }
    }

    print_time("performance_cpu.csv", "adi",
	       {"Ref", "Tiramisu"},
	       {median(duration_vector_1), median(duration_vector_2)});

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("adi", b_u_ref, b_u, 0.00001);

    if (PRINT_OUTPUT)
    {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_u);
        std::cout << "Reference " << std::endl;
        print_buffer(b_u_ref);
    }

    return 0;
}
