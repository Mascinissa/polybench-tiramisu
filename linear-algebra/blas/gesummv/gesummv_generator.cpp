#include <tiramisu/tiramisu.h>
#include "benchmarks.h"


using namespace tiramisu;

/*
   TODO description
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

    //Computations
    computation y_init("y_init", {i}, 0.0);
    computation y("y", {i,j}, p_float64);
    y.set_expression(y(i,j)+ A(i,j)*x(j)*alpha + B(i,j)*x(j)*beta);
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    y_init.then(y, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_input);
    buffer b_B("b_B", {N,N}, p_float64, a_input);
    buffer b_x("b_x", {N}, p_float64, a_input);
    buffer b_y("b_y", {N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    x.store_in(&b_x);
    

    //Store computations
    y_init.store_in(&b_y);
    y.store_in(&b_y,{i});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_x, &b_y}, "generated_gesummv.o");

    return 0;
}
