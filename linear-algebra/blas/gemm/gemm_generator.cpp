#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "gemm.h"


using namespace tiramisu;

/*
Generalized Matrix Multiply from BLAS.
It takes the following as inputs,
    • alpha, beta: scalars
    • A: PxQ matrix
    • B: QxR matrix
    • C: PxR matrix
and gives the following as output:
    • Cout: PxR array, where Cout = alpha*AB + beta*C
*/

int main(int argc, char **argv)
{
    tiramisu::init("gemm");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, P), j("j", 0, R), k("k", 0, Q);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {k, j}, p_float64);
    input C("C", {i, j}, p_float64);


    //Computations
    
    computation C_init("C_init", {i,j}, C(i,j)*beta);
    computation C_out("C_out", {i,k,j}, p_float64);
    C_out.set_expression(C(i,j)+A(i,k)*B(k,j)*alpha);
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    C_init.then(C_out, i);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {P,Q}, p_float64, a_input);
    buffer b_B("b_B", {Q,R}, p_float64, a_input);
    buffer b_C("b_C", {P,R}, p_float64, a_output);
    // buffer b_C_out("b_C_out", {P,R}, p_float64, a_output); 
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    

    //Store computations
    C_init.store_in(&b_C);
    C_out.store_in(&b_C, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C}, "generated_gemm.o");

    return 0;
}
