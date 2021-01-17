#include <tiramisu/tiramisu.h>
#include "benchmarks.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("trmm");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);
    constant MM("MM", M);

    //Iteration variables    
    var i("i", 0, M), j("j", 0, N), k("k", 0, M);
    

    //inputs
    input A("A", {i, k}, p_float64);
    input B("B", {i, j}, p_float64);


    //Computations
    

    computation AB("[MM,NN]->{AB[i,j,k]: 0<=i<MM and 0<=j<NN and i+1<=k<MM}", expr(), true, p_float64, global::get_implicit_function());
    AB.set_expression(AB(i,j,k) + A(k,i)*B(k,j));
    computation B_out("B_out", {i, j}, AB(i,j,0)*alpha);

    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    AB.then(B_out, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {M,M}, p_float64, a_input);
    buffer b_B("b_B", {M,N}, p_float64, a_output);
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    

    //Store computations
    AB.store_in(&b_B, {i,j});
    B_out.store_in(&b_B);

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B}, "generated_trmm.o");

    return 0;
}