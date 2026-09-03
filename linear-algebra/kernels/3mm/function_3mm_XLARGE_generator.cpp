#include <tiramisu/tiramisu.h>
#include <tiramisu/auto_scheduler/evaluator.h>
#include <tiramisu/auto_scheduler/search_method.h>
#include "function_3mm_XLARGE_wrapper.h"

using namespace tiramisu;

int main(int argc, char **argv)
{
    tiramisu::init("function_3mm_XLARGE");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, 1600), j("j", 0, 2400), k("k", 0, 2000), l("l", 0, 1800), m("m", 0, 2200);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {k, l}, p_float64);
    input C("C", {l, j}, p_float64);
    input D("D", {j, m}, p_float64);
    
    input AB_inp("AB_inp", {i, l}, p_float64);
    input CD_inp("CD_inp", {l, m}, p_float64);
    input E_inp("E_inp", {i,m}, p_float64);


    //Computations
    computation AB_init("AB_init", {i,l}, 0.0);
    computation AB("AB", {i,l,k}, AB_inp(i,l) + A(i,k)*B(k,l));

    computation CD_init("CD_init", {l,m}, 0.0);
    computation CD("CD", {l,m,j}, CD_inp(l,m) + C(l,j)*D(j,m));

    computation E_init("E_init", {i,m}, 0.0);
    computation E("E", {i,m,l}, E_inp(i,m) + AB_inp(i,l)*CD_inp(l,m));
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    AB_init.then(AB, l)
           .then(CD_init, computation::root)
           .then(CD, m)
           .then(E_init, computation::root)
           .then(E, m);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {1600,2000}, p_float64, a_input);
    buffer b_B("b_B", {2000,1800}, p_float64, a_input);
    buffer b_AB("b_AB", {1600,1800}, p_float64, a_temporary);
    buffer b_C("b_C", {1800,2400}, p_float64, a_input);
    buffer b_D("b_D", {2400,2200}, p_float64, a_input);
    buffer b_CD("b_CD", {1800,2200}, p_float64, a_temporary);
    buffer b_E("b_E", {1600,2200}, p_float64, a_output);
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    D.store_in(&b_D);
    AB_inp.store_in(&b_AB);
    CD_inp.store_in(&b_CD);
    E_inp.store_in(&b_E);
        

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
    tiramisu::codegen({&b_A, &b_B, &b_C, &b_D, &b_E}, "function_3mm_XLARGE.o");

    return 0;
}
