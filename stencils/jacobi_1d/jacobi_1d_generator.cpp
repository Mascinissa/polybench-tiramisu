#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "jacobi_1d.h"

using namespace tiramisu;

/*
Jacobi style stencil computation over 1D data with 3-point stencil pattern. 
*/

int main(int argc, char **argv)
{
    tiramisu::init("jacobi_1d");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i_f("i_f", 0, N);
    var t("t", 0, TSTEPS), i("i", 1, N-1);
    
    //inputs
    input A("A", {i_f}, p_float64);
    input B("B", {i_f}, p_float64);

    //Computations
    computation B_out("B_out", {t,i}, (A(i-1) + A(i) + A(i + 1))*0.33333);

    computation A_out("A_out", {t,i}, (B(i-1) + B(i) + B(i + 1))*0.33333);

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    B_out.then(A_out, t);
    
    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N}, p_float64, a_output);    
    buffer b_B("b_B", {N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);

    //Store computations
    A_out.store_in(&b_A, {i});
    B_out.store_in(&b_B, {i});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B}, "generated_jacobi_1d.o");

    return 0;
}
