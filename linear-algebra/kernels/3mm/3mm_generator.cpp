#include <tiramisu/tiramisu.h>
#include "benchmarks.h"
// #define eps 0.1

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("b_3mm");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, P), j("j", 0, S), k("k", 0, Q), l("l", 0, R), m("m", 0, T);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {k, l}, p_float64);
    input C("C", {l, j}, p_float64);
    input D("D", {j, m}, p_float64);


    //Computations
    computation AB_init("AB_init", {i,l}, 0.0);
    computation AB("AB", {i,l,k}, p_float64);
    AB.set_expression(AB(i,l,k) + A(i,k)*B(k,l));

    computation CD_init("CD_init", {l,m}, 0.0);
    computation CD("CD", {l,m,j}, p_float64);
    CD.set_expression(CD(l,m,j) + C(l,j)*D(j,m));

    computation E_init("E_init", {i,m}, 0.0);
    computation E("E", {i,m,l}, p_float64);
    E.set_expression(E(i,m,l) + AB(i,l,0)*CD(l,m,0));
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    AB_init.then(AB, computation::root)
           .then(CD_init, computation::root)
           .then(CD, computation::root)
           .then(E_init, computation::root)
           .then(E, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {P,Q}, p_float64, a_input);
    buffer b_B("b_B", {Q,R}, p_float64, a_input);
    buffer b_AB("b_AB", {P,R}, p_float64, a_temporary);
    buffer b_C("b_C", {R,S}, p_float64, a_input);
    buffer b_D("b_D", {S,T}, p_float64, a_input);
    buffer b_CD("b_CD", {R,T}, p_float64, a_temporary);
    buffer b_E("b_E", {P,T}, p_float64, a_output);
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    D.store_in(&b_D);
    

    //Store computations
    AB_init.store_in(&b_AB);
    CD_init.store_in(&b_CD);
    AB.store_in(&b_AB, {i,l});
    CD.store_in(&b_CD, {l,m});
    E_init.store_in(&b_E);
    E.store_in(&b_E, {i,m});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C, &b_D, &b_E}, "generated_3mm.o");

    return 0;
}
