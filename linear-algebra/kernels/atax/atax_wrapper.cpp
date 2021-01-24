#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_atax.o.h"
#include "polybench-tiramisu.h"
#include "atax.h"
#include <tiramisu/utils.h>

int atax_ref(Halide::Buffer<double> A, Halide::Buffer<double> x,
    Halide::Buffer<double> y)
{
    int i, j, k;
    Halide::Buffer<double> tmp(M);

    for (i = 0; i < N; i++)
        y(i) = 0;
    for (i = 0; i < M; i++) {
        tmp(i) = 0.0;
        for (j = 0; j < N; j++)
            tmp(i) = tmp(i) + A(i, j) * x(j);
        for (j = 0; j < N; j++)
            y(j) = y(j) + A(i, j) * tmp(i);
    }

    return 0;
}

int main(int argc, char** argv)
{
    std::vector<std::chrono::duration<double, std::milli> > duration_vector_1,
        duration_vector_2;
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
    Halide::Buffer<double> b_A(N, M), b_x(N), b_y(N), b_y_ref(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    // REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_array(b_A, b_x);

            transpose(b_A);
            auto start = std::chrono::high_resolution_clock::now();

            if (run_ref)
                atax_ref(b_A, b_x, b_y_ref);

            auto end = std::chrono::high_resolution_clock::now();
            transpose(b_A);

            duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_array(b_A, b_x);

            auto start = std::chrono::high_resolution_clock::now();
            if (run_tiramisu)
                atax(b_A.raw_buffer(), b_x.raw_buffer(), b_y.raw_buffer());

            auto end = std::chrono::high_resolution_clock::now();
            duration_vector_2.push_back(end - start);
        }
    }


    print_time("performance_cpu.csv", "atax", { "Ref", "Tiramisu" },
        { median(duration_vector_1), median(duration_vector_2) });

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("atax", b_y_ref, b_y, 0.0001);

    if (PRINT_OUTPUT) {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_y);
        std::cout << "Reference " << std::endl;
        print_buffer(b_y_ref);
    }

    return 0;
}
