#include <tiramisu/tiramisu.h>
#include "benchmarks.h"
// #define eps 0.1

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("b_2mm");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, P), j("j", 0, S), k("k", 0, Q), l("l", 0, R);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {k, l}, p_float64);
    input C("C", {l, j}, p_float64);
    input D("D", {i, j}, p_float64);


    //Computations
    computation D_init("D_init", {i,j}, D(i,j)*beta);
    computation D_out("D_out", {i,j,k,l}, p_float64);
    D_out.set_expression(D_out(i,j,k,l) + A(i,k)*B(k,l)*C(l,j)*alpha);
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    D_init.then(D_out, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {P,Q}, p_float64, a_input);
    buffer b_B("b_B", {Q,R}, p_float64, a_input);
    buffer b_C("b_C", {R,S}, p_float64, a_input);
    buffer b_D("b_D", {P,S}, p_float64, a_output);
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    D.store_in(&b_D);
    

    //Store computations
    D_init.store_in(&b_D);
    D_out.store_in(&b_D, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C, &b_D}, "generated_2mm.o");

    return 0;
}
