#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "syrk.h"

using namespace tiramisu;

/*
Symmetric rank k update.
It takes the following as inputs,
    • alpha, beta: scalars
    • A: NxM matrix
    • C: NxN symmetric matrix
and gives the following as output:
    • Cout: NxN matrix, where Cout = alpha*AAT + beta*C
*/

int main(int argc, char **argv)
{
    tiramisu::init("syrk");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);
    constant MM("MM", M);


    //Iteration variables    
    var i("i", 0, N), j("j"), k("k", 0, M);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input C("C", {i, j}, p_float64);


    //Computations
    computation C_beta("[NN]->{C_beta[i,j]: 0<=i<NN and 0<=j<=i}", expr(), true, p_float64, global::get_implicit_function());
    C_beta.set_expression(C(i,j)*beta);
    computation C_out("[MM,NN]->{C_out[i,k,j]: 0<=i<NN and 0<=j<=i and 0<=k<MM}", expr(), true, p_float64, global::get_implicit_function());
    C_out.set_expression(C(i,j)+ A(i,k)*A(j,k)*alpha);


    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    C_beta.then(C_out, i);


    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,M}, p_float64, a_input);
    buffer b_C("b_C", {N,N}, p_float64, a_output);


    //Store inputs
    A.store_in(&b_A);
    C.store_in(&b_C);


    //Store computations
    C_beta.store_in(&b_C);
    C_out.store_in(&b_C, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_C}, "generated_syrk.o");

    return 0;
}
