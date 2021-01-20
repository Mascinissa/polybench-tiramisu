#include <tiramisu/tiramisu.h>
#include "benchmarks.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("seidel_2d");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i_f("i_f", 0, N), j_f("j_f", 0, N);
    var t("t", 0, TSTEPS), i("i", 1, N-1), j("j", 1, N-1);
    
    //inputs
    input A("A", {i_f, j_f}, p_float64);

    //Computations
    computation A_out("A_out", {t,i,j}, (A(i-1, j-1) + A(i-1, j) + A(i-1, j+1)
		                               + A(i, j-1) + A(i, j) + A(i, j+1)
		                               + A(i+1, j-1) + A(i+1, j) + A(i+1, j+1))/9.0);

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    //no schedule

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_A("b_A", {N,N}, p_float64, a_output);    

    //Store inputs
    A.store_in(&b_A);


    //Store computations
    A_out.store_in(&b_A, {i,j});


    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_A}, "generated_seidel_2d.o");

    return 0;
}
