#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_mvt.o.h"
#include "polybench-tiramisu.h"
#include "mvt.h"
#include <tiramisu/utils.h>

int mvt_ref(Halide::Buffer<double> A, Halide::Buffer<double> y1, Halide::Buffer<double> y2, Halide::Buffer<double> x1, Halide::Buffer<double> x2)
{
    int i, j;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            x1(i) = x1(i) + A(i, j) * y1(j);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            x2(i) = x2(i) + A(j, i) * y2(j);

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
    Halide::Buffer<double> b_A(N, N), b_y1(N), b_y2(N), b_x1_ref(N), b_x2_ref(N), b_x1(N), b_x2(N);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    // REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_array(b_A, b_y1, b_y2 ,b_x1_ref, b_x2_ref);

            transpose(b_A);
            auto start = std::chrono::high_resolution_clock::now();

            if (run_ref)
                mvt_ref(b_A, b_y1, b_y2, b_x1_ref, b_x2_ref);

            auto end = std::chrono::high_resolution_clock::now();
            transpose(b_A);

            duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_array(b_A, b_y1, b_y2 ,b_x1, b_x2);


            auto start = std::chrono::high_resolution_clock::now();
            if (run_tiramisu)
                mvt(b_A.raw_buffer(), b_y1.raw_buffer(), b_y2.raw_buffer(), b_x1.raw_buffer(), b_x2.raw_buffer());

            auto end = std::chrono::high_resolution_clock::now();
            duration_vector_2.push_back(end - start);
        }
    }


    print_time("performance_cpu.csv", "mvt", { "Ref", "Tiramisu" },
        { median(duration_vector_1), median(duration_vector_2) });

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers_approximately("mvt x1", b_x1_ref, b_x1, 0.0001);
        compare_buffers_approximately("mvt x2", b_x2_ref, b_x2, 0.0001);

    if (PRINT_OUTPUT) {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_x1);
        std::cout << "Reference " << std::endl;
        print_buffer(b_x1_ref);
    }

    return 0;
}
