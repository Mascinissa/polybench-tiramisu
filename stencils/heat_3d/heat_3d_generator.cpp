#include <tiramisu/tiramisu.h>
#include "benchmarks.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("heat_3d");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i_f("i_f", 0, N), j_f("j_f", 0, N), k_f("k_f", 0, N);
    var t("t", 1, TSTEPS), i("i", 1, N-1), j("j", 1, N-1), k("k", 1, N-1);
    
    //inputs
    input A("A", {i_f, j_f, k_f}, p_float64);
    input B("B", {i_f, j_f, k_f}, p_float64);

    //Computations
    computation B_out("B_out", {t,i,j,k}, (A(i+1, j, k) - A(i, j, k)*2.0 + A(i-1, j, k))*0.125
                                        + (A(i, j+1, k) - A(i, j, k)*2.0 + A(i, j-1, k))*0.125
                                        + (A(i, j, k+1) - A(i, j, k)*2.0 + A(i, j, k-1))*0.125
                                        + A(i, j, k));

    computation A_out("A_out", {t,i,j,k}, (B(i+1, j, k) - B(i, j, k)*2.0 + B(i-1, j, k))*0.125
                                        + (B(i, j+1, k) - B(i, j, k)*2.0 + B(i, j-1, k))*0.125
                                        + (B(i, j, k+1) - B(i, j, k)*2.0 + B(i, j, k-1))*0.125
                                        + B(i, j, k));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    B_out.then(A_out, t);
    
    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N,N}, p_float64, a_output);    
    buffer b_B("b_B", {N,N,N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);

    //Store computations
    A_out.store_in(&b_A, {i,j,k});
    B_out.store_in(&b_B, {i,j,k});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B}, "generated_heat_3d.o");

    return 0;
}
