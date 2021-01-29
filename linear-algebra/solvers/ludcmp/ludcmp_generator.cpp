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
    var i("i", 0, N), j("j"), k("k"), l("l"), m("m"), n("n");
    

    //inputs
    input A("A", {i, i}, p_float64);
    input b("b", {i}, p_float64);


    //Computations
    computation A_sub("[NN]->{A_sub[i,j,k]: 0<=i<NN and 0<=j<i and 0<=k<j}", expr(), true, p_float64, global::get_implicit_function());
    A_sub.set_expression(A_sub(i,j,k) - A(i,k)*A(k,j));
    computation A_div("[NN]->{A_div[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    A_div.set_expression(A_sub(i,j,0)/A_sub(j,j,0));
    computation A_lu("[NN]->{A_lu[i,l,m]: 0<=i<NN and i<=l<NN and 0<=m<i}", expr(), true, p_float64, global::get_implicit_function());
    A_lu.set_expression(A_lu(i,l,m) - A_div(i,m)*A_div(m,l));

    computation y_init("y_init", {i}, b(i));
    computation y_temp("[NN]->{y_temp[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    y_temp.set_expression(y_temp(i,j) - A_lu(i, j, 0) * y_init(j));
    computation x_init("x_init", {i}, 0.0);
    computation y("[NN]->{y[i,n]: 0<=i<NN and NN-i<=n<NN}", expr(), true, p_float64, global::get_implicit_function());
    y.set_expression(y(i,n) - A(N-1-i, n) * x_init(n));
    computation x("[NN]->{x[i]: 0<=i<NN}", expr(), true, p_float64, global::get_implicit_function());
    x.set_expression(y(i,0) / A(N-1-i, N-1-i));
    

    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    A_sub.then(A_div,j)
         .then(A_lu, i)
         .then(y_init,computation::root)
         .then(y_temp,i)
         .then(x_init,computation::root)
         .then(y,computation::root)
         .then(x,i);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_output);    
    buffer b_b("b_b", {N}, p_float64, a_input);    
    buffer b_y("b_y", {N}, p_float64, a_output);    
    buffer b_x("b_x", {N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);    
    b.store_in(&b_b);    

    //Store computations
    A_sub.store_in(&b_A, {i,j});
    A_div.store_in(&b_A);
    A_lu.store_in(&b_A, {i,l});

    y_init.store_in(&b_y);
    y_temp.store_in(&b_y, {i});
    x_init.store_in(&b_x);
    y.set_access("[NN]->{y[i,n]->b_y[NN-1-i]}");
    x.set_access("[NN]->{x[i]->b_x[NN-1-i]}");


    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_b, &b_y, &b_x}, "generated_ludcmp.o");

    return 0;
}
