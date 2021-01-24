#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "correlation.h"
// #define eps 0.1

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    tiramisu::init("correlation");

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
    
    computation cov_init("cov_init", {i,j}, 0.0);
    computation cov("cov", {i,j,k}, p_float64);
    cov.set_expression(cov(i,j,k) + (data(k,i)-mean(0,i))*(data(k,j)-mean(0,j))/expr(cast(p_float64, N)));

    computation variance_init("variance_init", {i}, 0.0);
    computation variance("variance", {k, i}, p_float64);
    variance.set_expression(variance(k, i) + (data(k,i) - mean(0,i))*(data(k,i) - mean(0,i))/expr(cast(p_float64, N)));
    computation std("std", {i}, expr(o_sqrt, variance(0, i)));
    // computation std("std", {i}, expr(o_max, eps, expr(o_sqrt, variance(0, i))));

    computation corr("corr", {i, j}, cov(i,j,0)/(std(i)*std(j)));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    mean_init.then(mean, computation::root)
             .then(cov_init, computation::root)
             .then(cov, computation::root)
             .then(variance_init, computation::root)
             .then(variance, computation::root)
             .then(std, computation::root)
             .then(corr, computation::root);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_data("b_data", {N,M}, p_float64, a_input);
    buffer b_mean("b_mean", {M}, p_float64, a_temporary);
    buffer b_std("b_std", {M}, p_float64, a_temporary);
    buffer b_cov("b_cov", {M,M}, p_float64, a_temporary);   
    buffer b_corr("b_corr", {M,M}, p_float64, a_output);   
    

    //Store inputs
    data.store_in(&b_data);
    

    //Store computations
    mean_init.store_in(&b_mean);
    mean.store_in(&b_mean, {j});
    variance_init.store_in(&b_std);
    variance.store_in(&b_std, {i});
    std.store_in(&b_std);
    cov_init.store_in(&b_cov);
    cov.store_in(&b_cov, {i,j});
    corr.store_in(&b_corr);

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_data, &b_corr}, "generated_correlation.o");

    return 0;
}
