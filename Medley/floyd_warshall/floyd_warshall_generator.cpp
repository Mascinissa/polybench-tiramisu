#include <tiramisu/tiramisu.h>
#include "benchmarks.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("floyd_warshall");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, N), j("j", 0, N), k("k", 0, N);
    
    //inputs
    input paths("paths", {i, j}, p_float64);

    //Computations
    computation paths_update("paths_update", {k,i,j}, p_float64);
    paths_update.set_expression(expr(o_min, paths(i,j), paths(i,k) + paths(k,j)));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    // no_schedule
    
    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_paths("b_paths", {N,N}, p_float64, a_output);    

    //Store inputs
    paths.store_in(&b_paths);

    //Store computations
    paths_update.store_in(&b_paths, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_paths}, "generated_floyd_warshall.o");

    return 0;
}
