#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "trisolv.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("trisolv");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);

    //Iteration variables    
    var i("i", 0, N), j("j");
    

    //inputs
    input L("L", {i, i}, p_float64);
    input b("b", {i}, p_float64);


    //Computations

    computation x_init("x_init", {i}, b(i));
    computation x_sub("[NN]->{x_sub[i,j]: 0<=i<NN and 0<=j<i}", expr(), true, p_float64, global::get_implicit_function());
    x_sub.set_expression(x_sub(i,j) - L(i,j) * x_init(j));
    computation x_out("x_out", {i}, x_sub(i,0) / L(i,i));

    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    x_init.then(x_sub,i)
          .then(x_out, i);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_L("b_L", {N,N}, p_float64, a_input);    
    buffer b_b("b_b", {N}, p_float64, a_input);    
    buffer b_x("b_x", {N}, p_float64, a_output);    

    //Store inputs
    L.store_in(&b_L);    
    b.store_in(&b_b);    

    //Store computations
    x_init.store_in(&b_x);
    x_sub.store_in(&b_x, {i});
    x_out.store_in(&b_x);

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_L, &b_b, &b_x}, "generated_trisolv.o");

    return 0;
}
