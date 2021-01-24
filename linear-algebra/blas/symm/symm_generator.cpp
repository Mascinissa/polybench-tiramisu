#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "symm.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("symm");
    // function symm("symm");

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
    input C("C", {i, j}, p_float64);


    //Computations
    
    computation temp_init("temp_init", {i,j}, 0.0);
    computation temp("[MM,NN]->{temp[i,j,k]: 0<=i<MM and 0<=j<NN and 0<=k<i}", expr(), true, p_float64, global::get_implicit_function());
    temp.set_expression(temp(i,j,k)+A(i,k)*B(k,j));
    computation C_init("C_init", {i,j}, C(i, j)*beta+B(i, j)*A(i, i)*alpha+temp(i,j,0)*alpha);
    computation C_out("[MM,NN]->{C_out[i,k,j]: 0<=i<MM and 0<=k<i and 0<=j<NN}", expr(), true, p_float64, global::get_implicit_function());
    C_out.set_expression(C_out(i,k,j)+A(i,k)*B(k,j)*alpha);
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    temp_init.then(temp, computation::root)
             .then(C_init, computation::root)
             .then(C_out, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {M,M}, p_float64, a_input);
    buffer b_B("b_B", {M,N}, p_float64, a_input);
    buffer b_C("b_C", {M,N}, p_float64, a_output);
    buffer b_temp("b_temp", {M,N}, p_float64, a_temporary);
    // buffer b_C_out("b_C_out", {P,R}, p_float64, a_output); 
    

    //Store inputs
    A.store_in(&b_A);
    B.store_in(&b_B);
    C.store_in(&b_C);
    

    //Store computations
    temp_init.store_in(&b_temp);
    temp.store_in(&b_temp, {i,j});
    C_init.store_in(&b_C);
    C_out.store_in(&b_C, {k,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A, &b_B, &b_C}, "generated_symm.o");

    return 0;
}
