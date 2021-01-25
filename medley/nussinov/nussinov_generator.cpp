#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "nussinov.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("nussinov");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);

    //Iteration variables    
    var i("i", 0, N), j("j", 0, N), k("k");
    

    //inputs
    input table("table", {i, j}, p_float64);
    input seq("seq", {i}, p_float64);


    //Computations
    computation table_1("[NN]->{table_1[i,j]: 0<=i<NN and NN-i<=j<NN and 1<=j}", expr(), true, p_float64, global::get_implicit_function());
    table_1.set_expression(expr(o_max, table(N-1-i, j), table(N-1-i, j-1)));
    computation table_2("[NN]->{table_2[i,j]: 0<=i<NN and NN-i<=j<NN and 0<i}", expr(), true, p_float64, global::get_implicit_function());
    table_2.set_expression(expr(o_max, table(N-1-i, j), table(N-1-i+1, j)));
    computation table_3("[NN]->{table_3[i,j]: 0<=i<NN and NN-i<=j<NN and 1<=j and 0<i and NN-i<j}", expr(), true, p_float64, global::get_implicit_function());
    table_3.set_expression(expr(o_max, table(N-1-i, j), table(N-1-i+1, j-1)+((seq(N-1-i)+seq(j))==3.0)));
    computation table_4("[NN]->{table_4[i,j]: 0<=i<NN and NN-i<=j<NN and 1<=j and 0<i and j<=NN-i}", expr(), true, p_float64, global::get_implicit_function());
    table_4.set_expression(expr(o_max, table(N-1-i, j), table(N-1-i+1, j-1)));
    computation table_5("[NN]->{table_5[i,j,k]: 0<=i<NN and NN-i<=j<NN and NN-i<=k<j}", expr(), true, p_float64, global::get_implicit_function());
    table_5.set_expression(expr(o_max, table(N-1-i, j), table(N-1-i, k) + table(k+1, j)));
    
    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    table_1.then(table_2, j)
           .then(table_3, j)
           .then(table_4, j)
           .then(table_5, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_table("b_table", {N,N}, p_float64, a_output);    
    buffer b_seq("b_seq", {N}, p_float64, a_input);    

    //Store inputs
    table.store_in(&b_table);  
    seq.store_in(&b_seq);  

    //Store computations
    table_1.set_access("[NN]->{table_1[i,j]->b_table[NN-1-i,j]}");
    table_2.set_access("[NN]->{table_2[i,j]->b_table[NN-1-i,j]}");
    table_3.set_access("[NN]->{table_3[i,j]->b_table[NN-1-i,j]}");
    table_4.set_access("[NN]->{table_4[i,j]->b_table[NN-1-i,j]}");
    table_5.set_access("[NN]->{table_5[i,j,k]->b_table[NN-1-i,j]}");

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_table, &b_seq}, "generated_nussinov.o");

    return 0;
}
