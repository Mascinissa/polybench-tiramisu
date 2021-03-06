#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "cholesky.h"

using namespace tiramisu;

/*
Cholesky decomposition, which decomposes a matrix to triangular matrices. 
Only applicable when the input matrix is positive-definite.
It takes the following as input,
    • A: NxN positive-definite matrix
and gives the following as output:
    • L: NxN lower triangular matrix such that A = LL^T
*/

int main(int argc, char **argv)
{
    tiramisu::init("cholesky");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);

    //Iteration variables    
    var i("i"), j("j"), k("k"), l("l"), m("m");
    

    //inputs
    input A("A", {i, i}, p_float64);


    //Computations
    computation A_sub("[NN]->{A_sub[i,j,k]: 0<=i<NN and 0<=j<i and 0<=k<j}", expr(), true, p_float64, global::get_implicit_function());
    A_sub.set_expression(A_sub(i,j,k) - A(i,k)*A(j,k));
    computation A_div("[NN]->{A_div[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    A_div.set_expression(A_sub(i,j,0)/A_sub(j,j,0));
    computation A_diag("[NN]->{A_diag[i,l,m]: 0<=i<NN and l=i and 0<=m<i}", expr(), true, p_float64, global::get_implicit_function());
    A_diag.set_expression(A_diag(i,l,m) - A_div(i,m)*A_div(i,m));
    computation A_out("[NN]->{A_out[i,l]: 0<=i<NN and l=i}", expr(), true, p_float64, global::get_implicit_function());
    A_out.set_expression(expr(o_sqrt, A_diag(i,l,0)));

    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    A_sub.then(A_div,j)
         .then(A_diag, i)
         .then(A_out, i);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);    

    //Store computations
    A_sub.store_in(&b_A, {i,j});
    A_div.store_in(&b_A);
    A_diag.store_in(&b_A, {i,l});
    A_out.store_in(&b_A);

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A}, "generated_cholesky.o");

    return 0;
}
