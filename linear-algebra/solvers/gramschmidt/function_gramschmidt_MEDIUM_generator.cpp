#include <tiramisu/tiramisu.h>
#include <tiramisu/auto_scheduler/evaluator.h>
#include <tiramisu/auto_scheduler/search_method.h>
#include "function_gramschmidt_MEDIUM_wrapper.h"

using namespace tiramisu;


int main(int argc, char **argv)
{
    tiramisu::init("function_gramschmidt_MEDIUM");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, 200), j("j", 0, 240), k("k", 0, 240), l("l"), m("m");
    
    //inputs
    input A("A", {i, k}, p_float64);
    input Q("Q", {i, k}, p_float64);
    input R("R", {j, j}, p_float64);
    
    //Computations
    // The column norm is accumulated directly into R(k,k) (one distinct
    // location per k) instead of a single shared scalar buffer: with a
    // shared scalar, Tiramisu generated code whose per-k reset did not
    // take effect (nrm accumulated across all k), diverging from
    // PolyBench. See README "Validation status".
    computation R_diag_init("{R_diag_init[k]: 0<=k<240}", expr(), true, p_float64, global::get_implicit_function());
    R_diag_init.set_expression(0.0);

    computation R_diag_acc("{R_diag_acc[k,i]: 0<=k<240 and 0<=i<200}", expr(), true, p_float64, global::get_implicit_function());
    R_diag_acc.set_expression(R(k,k) + A(i, k) * A(i, k));

    computation R_diag("{R_diag[k]: 0<=k<240}", expr(), true, p_float64, global::get_implicit_function());
    R_diag.set_expression(expr(o_sqrt, R(k,k)));

    computation Q_out("{Q_out[k,i]: 0<=k<240 and 0<=i<200}", expr(), true, p_float64, global::get_implicit_function());
    Q_out.set_expression(A(i,k) / R(k,k));

    computation R_up_init("{R_up_init[k,j]: 0<=k<240 and k+1<=j<240}", expr(), true, p_float64, global::get_implicit_function());
    R_up_init.set_expression(0.0);

    computation R_up("{R_up[k,j,i]: 0<=k<240 and k+1<=j<240 and 0<=i<200}", expr(), true, p_float64, global::get_implicit_function());
    R_up.set_expression(R(k,j) + Q(i,k) * A(i, j)); 

    computation A_out("{A_out[k,j,i]: 0<=k<240 and k+1<=j<240 and 0<=i<200}", expr(), true, p_float64, global::get_implicit_function());
    A_out.set_expression(A(i,j) - Q(i,k) * R(k, j));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    R_diag_init.then(R_diag_acc, k)
            .then(R_diag, k)
            .then(Q_out, k)
            .then(R_up_init, k)
            .then(R_up, j)
            .then(A_out, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {200,240}, p_float64, a_output);
    buffer b_R("b_R", {240,240}, p_float64, a_output);
    buffer b_Q("b_Q", {200,240}, p_float64, a_output);  

    //Store inputs
    A.store_in(&b_A);    
    Q.store_in(&b_Q);    
    R.store_in(&b_R);    
    
    //Store computations
    R_diag_init.store_in(&b_R, {k,k});
    R_diag_acc.store_in(&b_R, {k,k});
    R_diag.store_in(&b_R, {k,k});
    Q_out.store_in(&b_Q, {i,k});
    R_up_init.store_in(&b_R);
    R_up.store_in(&b_R, {k,j});
    A_out.store_in(&b_A, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_Q, &b_R}, "function_gramschmidt_MEDIUM.o");

    return 0;
}
