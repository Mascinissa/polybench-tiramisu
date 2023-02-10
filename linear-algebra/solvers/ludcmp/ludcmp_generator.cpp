#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "ludcmp.h"

using namespace tiramisu;

/*
This kernel solves a system of linear equations using LU decomposition followed by forward and backward
substitutions.
It takes the following as inputs,
    • A: NxN matrix
    • b: vector of length N
and gives the following as output:
    • x: vector of length N, where Ax = b
The matrix A is first decomposed into L and U using the same algorithm as in lu. Then the two
triangular systems are solved to find x.
*/

int main(int argc, char **argv)
{
    tiramisu::init("ludcmp");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);

    //Iteration variables    
    var i("i", 0, N), j("j"), k("k"), l("l"), m("m"), n("n"), i_rev("i_rev");
    

    //inputs
    input A("A", {i, i}, p_float64);
    input b("b", {i}, p_float64);
    input w("w", {}, p_float64);
    input y("y", {i}, p_float64);
    input x("x", {i}, p_float64);


    //Computations
    computation w_init("[NN]->{w_init[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    w_init.set_expression(A(i,j));

    computation w_sub("[NN]->{w_sub[i,j,k]: 0<=i<NN and 0<=j<i and 0<=k<j}", expr(), true, p_float64, global::get_implicit_function());
    w_sub.set_expression(w(0) - A(i,k)*A(k,j));

    computation A_div("[NN]->{A_div[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    A_div.set_expression(w(0)/A(j,j));

    computation w_init2("[NN]->{w_init2[i,j]: 0<=i<NN and i<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    w_init2.set_expression(A(i,j));

    computation w_sub2("[NN]->{w_sub2[i,j,k]: 0<=i<NN and i<=j<NN and 0<=k<i}", expr(), true, p_float64, global::get_implicit_function());
    w_sub2.set_expression(w(0) - A(i,k)*A(k,j));

    computation A_cpy("[NN]->{A_cpy[i,j]: 0<=i<NN and i<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    A_cpy.set_expression(w(0));

    computation w_init3("[NN]->{w_init3[i]: 0<=i<NN}", expr(), true, p_float64, global::get_implicit_function());
    w_init3.set_expression(b(i));

    computation w_sub3("[NN]->{w_sub3[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    w_sub3.set_expression(w(0) - A(i,j)*y(j));

    computation y_init("[NN]->{y_init[i]: 0<=i<NN}", expr(), true, p_float64, global::get_implicit_function());
    y_init.set_expression(w(0));

    computation w_init4("[NN]->{w_init4[i]: 0<=i<NN}", expr(), true, p_float64, global::get_implicit_function());
    w_init4.set_expression(y(i));

    computation w_sub4("[NN]->{w_sub4[i,j]: 0<=i<NN and i+1<=j<NN }", expr(), true, p_float64, global::get_implicit_function());
    w_sub4.set_expression(w(0) - A(i,j)*x(j));

    computation x_comp("[NN]->{x_comp[i]: 0<=i<NN}", expr(), true, p_float64, global::get_implicit_function());
    x_comp.set_expression(w(0)/A(i,i));

    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    w_init4.loop_reversal(i, i_rev);
    w_sub4.loop_reversal(i, i_rev);
    x_comp.loop_reversal(i, i_rev);

    w_init.then(w_sub, j)
          .then(A_div, j)
          .then(w_init2, i)
          .then(w_sub2, j)
          .then(A_cpy, j)
          .then(w_init3, computation::root)
          .then(w_sub3, i)
          .then(y_init, i)
          .then(w_init4, computation::root)
          .then(w_sub4, i_rev)
          .then(x_comp, i_rev);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_output);    
    buffer b_b("b_b", {N}, p_float64, a_input);    
    buffer b_y("b_y", {N}, p_float64, a_output);
    buffer b_x("b_x", {N}, p_float64, a_output);
    buffer b_w("b_w", {1}, p_float64, a_temporary);

    //Store inputs
    A.store_in(&b_A);    
    b.store_in(&b_b);
    w.store_in(&b_w);
    y.store_in(&b_y);
    x.store_in(&b_x);

    //Store computations
    w_init.store_in(&b_w, {});
    w_init2.store_in(&b_w, {});
    w_init3.store_in(&b_w, {});
    w_init4.store_in(&b_w, {});
    w_sub.store_in(&b_w, {});
    w_sub2.store_in(&b_w, {});
    w_sub3.store_in(&b_w, {});
    w_sub4.store_in(&b_w, {});
    A_div.store_in(&b_A);
    A_cpy.store_in(&b_A);

    y_init.store_in(&b_y);
    x_comp.store_in(&b_x);



    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_b, &b_y, &b_x}, "generated_ludcmp.o");

    return 0;
}
