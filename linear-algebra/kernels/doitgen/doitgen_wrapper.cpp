#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_doitgen.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>

int doitgen_ref(Halide::Buffer<double> A, Halide::Buffer<double> x)
{
int r, q, p, s;
Halide::Buffer<double> sum(P);

for (r = 0; r < R; r++)
  for (q = 0; q < Q; q++) {
    for (p = 0; p < P; p++) {
      sum(p) = 0.0;
      for (s = 0; s < P; s++)
        sum(p) += A(r, q, s) * x(s, p);
    }
    for (p = 0; p < P; p++)
      A(r, q, p) = sum(p);
  }

    return 0;
}

int main(int argc, char** argv)
{
    std::vector<std::chrono::duration<double, std::milli> > duration_vector_1, duration_vector_2;
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
    Halide::Buffer<double> b_A(P, Q, R), b_A_ref(P, Q, R), b_x(P, P);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    // REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_buffer(b_A_ref, (double)2);
            init_buffer(b_x, (double)4);

            transpose(b_A_ref);
            transpose(b_x);
            auto start = std::chrono::high_resolution_clock::now();

            if (run_ref)
                doitgen_ref(b_A_ref, b_x);

            auto end = std::chrono::high_resolution_clock::now();
            transpose(b_A_ref);
            transpose(b_x);

            duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_buffer(b_A, (double)2);
            init_buffer(b_x, (double)4);

            auto start = std::chrono::high_resolution_clock::now();
            if (run_tiramisu)
                doitgen(b_A.raw_buffer(), b_x.raw_buffer());

            auto end = std::chrono::high_resolution_clock::now();
            duration_vector_2.push_back(end - start);
        }
    }


    print_time("performance_cpu.csv", "doitgen", { "Ref", "Tiramisu" },
        { median(duration_vector_1), median(duration_vector_2) });

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("doitgen", b_A_ref, b_A);

    if (PRINT_OUTPUT) {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_A);
        std::cout << "Reference " << std::endl;
        print_buffer(b_A_ref);
    }

    return 0;
}
