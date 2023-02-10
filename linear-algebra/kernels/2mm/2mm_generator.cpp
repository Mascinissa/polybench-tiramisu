#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "2mm.h"

using namespace tiramisu;

/*
Linear algebra kernel that consists of two matrix multiplications.
It takes the following as inputs,
    • alpha, beta: scalars
    • A: PxQ matrix
    • B: QxR matrix
    • C: RxS matrix
    • D: PxS matrix
and gives the following as output:
    • E: PxS matrix, where E = alpha*ABC + beta*D
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
    input tmp("tmp", {i,l}, p_float64);

    //Computations
    computation tmp_init("tmp_init",{i,l}, 0.0);
    computation tmp_prod("tmp_prod",{i,l,k}, tmp(i,l) + A(i,k)*B(k,l)*alpha);

    computation D_beta("D_beta", {i,j}, D(i,j)*beta);
    computation D_prod("D_prod", {i,j,l}, D(i,j)+tmp(i,l)*C(l,j));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    tmp_init.then(tmp_prod,l)
            .then(D_beta, computation::root)
            .then(D_prod, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {P,Q}, p_float64, a_input);
    buffer b_B("b_B", {Q,R}, p_float64, a_input);
    buffer b_C("b_C", {R,S}, p_float64, a_input);
    buffer b_D("b_D", {P,S}, p_float64, a_output);
    buffer b_tmp("b_tmp", {P,R}, p_float64, a_temporary);


    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    D.store_in(&b_D);
    tmp.store_in(&b_tmp);


    //Store computations
    tmp_init.store_in(&b_tmp);
    tmp_prod.store_in(&b_tmp, {i,l});
    D_beta.store_in(&b_D);
    D_prod.store_in(&b_D, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C, &b_D}, "generated_2mm.o");

    return 0;
}
