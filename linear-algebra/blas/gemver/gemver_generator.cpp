#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "gemver.h"

using namespace tiramisu;

/*
Multiple matrix-vector multiplication.
It takes the following as inputs,
    • alpha, beta: scalars
    • A: NxN matrix
    • u1, u2, v1, v2, y, z: vectors of length N
and gives the following as outputs:
    • A’: NxN matrix, where A' = A + u1.v1 + u2.v2
    • x: vector of length N, where x = beta*A'*y + z
    • w: vector of length N, where w = alpha*A'*x
*/

int main(int argc, char **argv)
{
    tiramisu::init("gemver");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, N), j("j", 0, N);
    

    //inputs
    input A("A", {i, j}, p_float64);
    input u1("u1", {i}, p_float64);
    input u2("u2", {i}, p_float64);
    input v1("v1", {i}, p_float64);
    input v2("v2", {i}, p_float64);
    input y("y", {i}, p_float64);
    input z("z", {i}, p_float64);


    //Computations
    
    computation A_hat("A_hat", {i,j}, A(i, j) + u1(i)*v1(j) + u2(i)*v2(j));
    computation x_temp("x_temp", {i,j}, p_float64);
    x_temp.set_expression(x_temp(i,j) + A_hat(j, i)*y(j)*beta);
    computation x("x", {i}, x_temp(i, 0) + z(i));
    computation w("w", {i,j}, p_float64);
    w.set_expression(w(i,j) + A_hat(i, j) * x(j)*alpha);


    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    A_hat.then(x_temp, computation::root)
         .then(x, computation::root)
         .then(w, computation::root);


    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_input);
    buffer b_u1("b_u1", {N}, p_float64, a_input);
    buffer b_u2("b_u2", {N}, p_float64, a_input);
    buffer b_v1("b_v1", {N}, p_float64, a_input);
    buffer b_v2("b_v2", {N}, p_float64, a_input);
    buffer b_z("b_z", {N}, p_float64, a_input);
    buffer b_y("b_y", {N}, p_float64, a_input);
    buffer b_A_hat("b_A_hat", {N,N}, p_float64, a_output);
    buffer b_x("b_x", {N}, p_float64, a_output);
    buffer b_w("b_w", {N}, p_float64, a_output);

    //Store inputs
    A.store_in(&b_A);
    u1.store_in(&b_u1);
    u2.store_in(&b_u2);
    v1.store_in(&b_v1);
    v2.store_in(&b_v2);
    y.store_in(&b_y);
    z.store_in(&b_z);
    
    //Store computations
    A_hat.store_in(&b_A_hat);
    x_temp.store_in(&b_x, {i});
    x.store_in(&b_x);
    w.store_in(&b_w, {i});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_u1, &b_u2, &b_v1, &b_v2, &b_y, &b_z, &b_A_hat, &b_x, &b_w}, "generated_gemver.o");

    return 0;
}
