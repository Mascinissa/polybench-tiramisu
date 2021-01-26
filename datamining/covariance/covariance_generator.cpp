#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "covariance.h"

using namespace tiramisu;

/*
Computes the covariance, a measure from statistics that show how linearly related two variables are.
It takes the following as input,
    •data:NxMmatrix that representsNdata points, each with M attributes,
and gives the following as output:
    •cov:MxMmatrix where the i, j-th element is the covariance between i and j.  The matrix issymmetric.          
*/

int main(int argc, char **argv)
{
    tiramisu::init("covariance");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 0, M), j("j", 0, M), k("k", 0, N), l("l", 0, N);
    

    //inputs
    input data("data", {l, j}, p_float64);


    //Computations
    
    computation mean_init("mean_init", {j}, 0.0);
    computation mean("mean", {l,j}, p_float64);
    mean.set_expression(mean(l,j) + data(l,j)/expr(cast(p_float64, N)));
    
    computation cov_init("conv_init", {i,j}, 0.0);
    computation cov("cov", {i,j,k}, p_float64);
    cov.set_expression(cov(i,j,k) + (data(k,i)-mean(0,i))*(data(k,j)-mean(0,j))/expr(cast(p_float64, N-1)));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    mean_init.then(mean, computation::root)
             .then(cov_init, computation::root)
             .then(cov, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_data("b_data", {N,M}, p_float64, a_input);
    buffer b_mean("b_mean", {M}, p_float64, a_temporary);
    buffer b_cov("b_cov", {M,M}, p_float64, a_output);   
    

    //Store inputs
    data.store_in(&b_data);
    

    //Store computations
    mean_init.store_in(&b_mean);
    mean.store_in(&b_mean, {j});
    cov_init.store_in(&b_cov);
    cov.store_in(&b_cov, {i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_data, &b_cov}, "generated_covariance.o");

    return 0;
}
