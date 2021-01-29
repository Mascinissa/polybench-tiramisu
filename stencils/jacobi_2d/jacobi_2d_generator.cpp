#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "jacobi_2d.h"

using namespace tiramisu;

/*
Jacobi-style stencil computation over 2D data with 5-point stencil pattern. The computation is simplified
as simply taking the average of five points.
*/

int main(int argc, char **argv)
{
    tiramisu::init("jacobi_2d");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i_f("i_f", 0, N), j_f("j_f", 0, N);
    var t("t", 0, TSTEPS), i("i", 1, N-1), j("j", 1, N-1);
    
    //inputs
    input A("A", {i_f, j_f}, p_float64);
    input B("B", {i_f, j_f}, p_float64);

    //Computations
    computation B_out("B_out", {t,i,j}, (A(i, j) + A(i, j-1) + A(i, 1+j) + A(1+i, j) + A(i-1, j))*0.2);

    computation A_out("A_out", {t,i,j}, (B(i, j) + B(i, j-1) + B(i, 1+j) + B(1+i, j) + B(i-1, j))*0.2);

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    B_out.then(A_out, t);
    
    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_output);    
    buffer b_B("b_B", {N,N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);

    //Store computations
    A_out.store_in(&b_A, {i,j});
    B_out.store_in(&b_B, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B}, "generated_jacobi_2d.o");

    return 0;
}
