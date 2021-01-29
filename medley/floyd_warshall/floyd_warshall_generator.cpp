#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "floyd_warshall.h"


using namespace tiramisu;

/*
Floyd-Warshall computes the shortest paths between each pair of nodes in a graph.
It takes the following as input,
    • paths: NxN matrix, where the i, jth entry represents the cost of taking an edge from i to j. Set to infinity
if there is no edge connecting i to j.
and gives the following as output:
    • paths: NxN matrix, where the i, jth entry represents the shortest path length from i to j
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
