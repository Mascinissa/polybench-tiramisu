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
    input mean("mean", {j}, p_float64);
    input cov("cov", {i,j}, p_float64);

    //Computations
    
    computation mean_init("mean_init", {j}, 0.0);
    computation mean_sum("mean_sum", {j,l}, mean(j) + data(l,j));

    computation mean_div("mean_div", {j}, mean(j) /expr(cast(p_float64, N)));

    computation data_sub("data_sub", {l,j}, data(l,j)-mean(j));

    computation cov_init("conv_init", {i,j}, 0.0);

    computation cov_prod("cov_prod", {i,j,k}, cov(i,j) + data(k,i)*data(k,j));

    computation cov_div("cov_div", {i,j}, cov(i,j)/expr(cast(p_float64, N-1)) );

    computation cov_sym("cov_sym", {i,j}, cov(i,j));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    mean_init.then(mean_sum, j)
             .then(mean_div,j)
             .then(data_sub,computation::root)
             .then(cov_init, computation::root)
             .then(cov_prod, j)
             .then(cov_div, j)
             .then(cov_sym, j);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_data("b_data", {N,M}, p_float64, a_input);
    buffer b_mean("b_mean", {M}, p_float64, a_temporary);
    buffer b_cov("b_cov", {M,M}, p_float64, a_output);   
    

    //Store inputs
    data.store_in(&b_data);
    mean.store_in(&b_mean);
    cov.store_in(&b_cov);


    //Store computations
    mean_init.store_in(&b_mean);
    mean_sum.store_in(&b_mean, {j});
    mean_div.store_in(&b_mean, {j});
    data_sub.store_in(&b_data);
    cov_init.store_in(&b_cov);
    cov_prod.store_in(&b_cov, {i,j});
    cov_div.store_in(&b_cov, {i,j});
    cov_sym.store_in(&b_cov, {j,i});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_data, &b_cov}, "generated_covariance.o");

    return 0;
}
