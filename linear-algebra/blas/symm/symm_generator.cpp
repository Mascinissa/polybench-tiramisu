#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "symm.h"

using namespace tiramisu;

/*
Symmetric matrix matrix multiplication from BLAS.
It takes the following as inputs,
    • alpha, beta: scalars
    • A: M x M symmetric matrix
    • B,C: M x N matrices
and gives the following as output:
    • Cout: M x N matrix, where Cout = alpha*AB + beta
*/

int main(int argc, char **argv)
{
    tiramisu::init("symm");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);
    constant MM("MM", M);

    //Iteration variables    
    var i("i", 0, M), j("j", 0, N), k("k", 0, M);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {i, j}, p_float64);
    input C("C", {i, j}, p_float64);
    input temp("temp", {}, p_float64);


    //Computations
//    computation temp_init("temp_init", {i,j}, 0.0);
    computation temp_init("[MM,NN]->{temp_init[i,j]: 0<=i<MM and 0<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    temp_init.set_expression(0.0);
    computation temp_comp("[MM,NN]->{temp_comp[i,j,k]: 0<=i<MM and 0<=j<NN and 0<=k<i}", expr(), true, p_float64, global::get_implicit_function());
    temp_comp.set_expression(temp(0)+A(i,k)*B(k,j));
    computation C_r("[MM,NN]->{C_r[i,j,k]: 0<=i<MM and 0<=j<NN and 0<=k<i}", expr(), true, p_float64, global::get_implicit_function());
    C_r.set_expression(C(k,j) + B(i,j)*A(i,k)*alpha);
    computation C_out("[MM,NN]->{C_out[i,j]: 0<=i<MM and 0<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    C_out.set_expression(C(i, j)*beta + B(i, j)*A(i, i)*alpha +  temp(0)*alpha);
    
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    temp_init.then(C_r, j)
             .then(temp_comp, k)
             .then(C_out, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {M,M}, p_float64, a_input);
    buffer b_B("b_B", {M,N}, p_float64, a_input);
    buffer b_C("b_C", {M,N}, p_float64, a_output);
    buffer b_temp("b_temp", {1}, p_float64, a_temporary);

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    temp.store_in(&b_temp);
    

    //Store computations
    temp_init.store_in(&b_temp, {});
    temp_comp.store_in(&b_temp, {});
    C_r.store_in(&b_C, {k,j});
    C_out.store_in(&b_C);

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C}, "generated_symm.o");

    return 0;
}
