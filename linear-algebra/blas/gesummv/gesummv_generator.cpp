#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "gesummv.h"


using namespace tiramisu;

/*
Summed matrix-vector multiplications.
It takes the following as inputs,
    • alpha, beta: scalars
    • A, B: N x N matrix
    • x: vector of length N
and gives the following as outputs:
    • y: vector of length N, where y = alpha*Ax + beta*Bx
*/

int main(int argc, char **argv)
{
    tiramisu::init("gesummv");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, N), j("j", 0, N);
    

    //inputs
    input A("A", {i, j}, p_float64);
    input B("B", {i, j}, p_float64);
    input x("x", {i}, p_float64);
    input y("y", {i}, p_float64);
    input tmp("tmp", {i}, p_float64);

    //Computations
    computation tmp_init("tmp_init", {i}, 0.0);
    computation y_init("y_init", {i}, 0.0);
    computation tmp_comp("tmp_comp", {i,j}, tmp(i)+A(i,j)*x(j));
    computation y_comp1("y_comp1", {i,j}, y(i)+B(i,j)*x(j));
    computation y_comp2("y_comp2", {i}, tmp(i)*alpha + y(i)*beta);
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    tmp_init.then(y_init, i)
            .then(tmp_comp,i)
            .then(y_comp1,{j})
            .then(y_comp2,i);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_tmp("b_tmp", {N}, p_float64, a_temporary);
    buffer b_A("b_A", {N,N}, p_float64, a_input);
    buffer b_B("b_B", {N,N}, p_float64, a_input);
    buffer b_x("b_x", {N}, p_float64, a_input);
    buffer b_y("b_y", {N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    x.store_in(&b_x);
    y.store_in(&b_y);
    tmp.store_in(&b_tmp);

    //Store computations
    tmp_init.store_in(&b_tmp);
    tmp_comp.store_in(&b_tmp,{i});
    y_init.store_in(&b_y);
    y_comp1.store_in(&b_y,{i});
    y_comp2.store_in(&b_y);


    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_x, &b_y}, "generated_gesummv.o");

    return 0;
}
