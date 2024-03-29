#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "gramschmidt.h"

using namespace tiramisu;

/*
QR Decomposition with Modified Gram Schmidt.
It takes the following as input,
    • A: MxN rank N matrix (M>N).
and gives the following as outputs:
    • Q: MxN orthogonal matrix
    • R: NxN upper triangular matrix
such that A = QR.
*/

int main(int argc, char **argv)
{
    tiramisu::init("gramschmidt");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N), MM("MM", M);

    //Iteration variables    
    var i("i", 0, M), j("j", 0, N), k("k", 0, N), l("l"), m("m");
    
    //inputs
    input A("A", {i, k}, p_float64);
    input Q("Q", {i, k}, p_float64);
    input R("R", {j, j}, p_float64);
    input nrm("nrm", {}, p_float64);

    //Computations
    computation nrm_init("[NN]->{nrm_init[k]: 0<=k<NN}", expr(), true, p_float64, global::get_implicit_function());
    nrm_init.set_expression(eps);

    computation nrm_comp("[NN,MM]->{nrm_comp[k,i]: 0<=k<NN and 0<=i<MM}", expr(), true, p_float64, global::get_implicit_function());
    nrm_comp.set_expression(nrm(0) + A(i, k) * A(i, k));

    computation R_diag("[NN]->{R_diag[k]: 0<=k<NN}", expr(), true, p_float64, global::get_implicit_function());
    R_diag.set_expression(expr(o_sqrt, nrm(0)));

    computation Q_out("[NN,MM]->{Q_out[k,i]: 0<=k<NN and 0<=i<MM}", expr(), true, p_float64, global::get_implicit_function());
    Q_out.set_expression(A(i,k) / R(k,k));

    computation R_up_init("[NN]->{R_up_init[k,j]: 0<=k<NN and k+1<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    R_up_init.set_expression(0.0);

    computation R_up("[NN,MM]->{R_up[k,j,i]: 0<=k<NN and k+1<=j<NN and 0<=i<MM}", expr(), true, p_float64, global::get_implicit_function());
    R_up.set_expression(R(k,j) + Q(i,k) * A(i, j)); 

    computation A_out("[NN,MM]->{A_out[k,j,i]: 0<=k<NN and k+1<=j<NN and 0<=i<MM}", expr(), true, p_float64, global::get_implicit_function());
    A_out.set_expression(A(i,j) - Q(i,k) * R(k, j)) ;

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    nrm_init.then(nrm_comp, k)
            .then(R_diag, k)
            .then(Q_out, k)
            .then(R_up_init, k)
            .then(R_up, j)
            .then(A_out, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {M,N}, p_float64, a_output);
    buffer b_nrm("b_nrm", {1}, p_float64, a_temporary);
    buffer b_R("b_R", {N,N}, p_float64, a_output);
    buffer b_Q("b_Q", {M,N}, p_float64, a_output);  

    //Store inputs
    A.store_in(&b_A);    
    Q.store_in(&b_Q);    
    R.store_in(&b_R);
    nrm.store_in(&b_nrm);

    //Store computations
    nrm_init.store_in(&b_nrm, {});
    nrm_comp.store_in(&b_nrm, {});
    R_diag.store_in(&b_R, {k,k});
    Q_out.store_in(&b_Q, {i,k});
    R_up_init.store_in(&b_R);
    R_up.store_in(&b_R, {k,j});
    A_out.store_in(&b_A, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_Q, &b_R}, "generated_gramschmidt.o");

    return 0;
}
