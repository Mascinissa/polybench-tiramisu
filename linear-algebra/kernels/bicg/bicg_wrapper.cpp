#include <Halide.h>
#include <tiramisu/tiramisu.h>
#include <iostream>
#include "generated_bicg.o.h"
#include "benchmarks.h"
#include <tiramisu/utils.h>

int bicg_ref(Halide::Buffer<double> A, Halide::Buffer<double> p, Halide::Buffer<double> r, Halide::Buffer<double> q, Halide::Buffer<double> s)
{
    int i, j, k;

    for (i = 0; i < M; i++)
        s(i) = 0;
    for (i = 0; i < N; i++){
        q(i) = 0.0;
        for (j = 0; j < M; j++){
            s(j) = s(j) + r(i) * A(i, j);
            q(i) = q(i) + A(i, j) * p(j);
        }
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
    Halide::Buffer<double> b_A(M, N), b_p(M), b_r(N), b_q_ref(N), b_s_ref(M), b_q(N), b_s(M);
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    // REFERENCE
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_buffer(b_q_ref, (double)1);
            init_buffer(b_s_ref, (double)2);
            init_buffer(b_A, (double)3);
            init_buffer(b_p, (double)4);
            init_buffer(b_r, (double)5);

            transpose(b_A);
            auto start = std::chrono::high_resolution_clock::now();

            if (run_ref)
                bicg_ref(b_A, b_p, b_r, b_q_ref, b_s_ref);

            auto end = std::chrono::high_resolution_clock::now();
            transpose(b_A);

            duration_vector_1.push_back(end - start);
        }
    }

    // TIRAMISU
    {
        for (int i = 0; i < NB_TESTS; ++i) {
            init_buffer(b_q, (double)1);
            init_buffer(b_s, (double)2);
            init_buffer(b_A, (double)3);
            init_buffer(b_p, (double)4);
            init_buffer(b_r, (double)5);

            auto start = std::chrono::high_resolution_clock::now();
            if (run_tiramisu)
                bicg(b_A.raw_buffer(), b_p.raw_buffer(), b_r.raw_buffer(), b_q.raw_buffer(), b_s.raw_buffer());

            auto end = std::chrono::high_resolution_clock::now();
            duration_vector_2.push_back(end - start);
        }
    }


    print_time("performance_cpu.csv", "bicg", { "Ref", "Tiramisu" },
        { median(duration_vector_1), median(duration_vector_2) });

    if (CHECK_CORRECTNESS && run_ref && run_tiramisu)
        compare_buffers("bicg q", b_q_ref, b_q);
        compare_buffers("bicg s", b_s_ref, b_s);

    if (PRINT_OUTPUT) {
        std::cout << "Tiramisu " << std::endl;
        print_buffer(b_q);
        std::cout << "Reference " << std::endl;
        print_buffer(b_q_ref);
    }

    return 0;
}
