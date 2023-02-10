#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "durbin.h"

using namespace tiramisu;

/*
Durbing is an algorithm for solving Yule-Walker equations, which is a special case of Toelitz systems.
It takes the following as input,
   • r: vector of length N.
and gives the following as output:
   • y: vector of length N
such that Ty = -r where T is a symmetric, unit-diagonal, Toeplitz matrix defined by the vector [1, r_0, ..., r_(N-1)].
*/

int main(int argc, char **argv)
{
    tiramisu::init("durbin");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 
    constant NN("NN", N);

    //Iteration variables    
    var i("i", 0, N), j("j", 0, N), k("k", 1, N);
    var d("d",0, 1);
    
    //inputs
    input r("r", {i}, p_float64);
    input y("y", {i}, p_float64);
    input alpha("alpha", {}, p_float64);
    input beta("beta", {}, p_float64);
    input sum("sum", {}, p_float64);

    computation y_0("y_0", {d}, -r(0));
    computation beta_0("beta_0", {d}, 1.0);
    computation alpha_0("alpha_0", {d}, -r(0));

    //Computations
    computation beta_1("beta_1", {k}, (1-alpha(0)*alpha(0))*beta(0));
    computation sum_1("sum_1", {k}, 0.0);
    computation sum_2("[NN]->{sum_2[k,i]: 1<=k<NN and 0<=i<k}", expr(), true, p_float64, global::get_implicit_function());
    sum_2.set_expression(sum(0)+r(k-i-1)*y(i));
    computation alpha_1("alpha_1", {k}, - (r(k) + sum(0))/beta(0));
    computation z("[NN]->{z[k,i]: 1<=k<NN and 0<=i<k}", expr(), true, p_float64, global::get_implicit_function());
    z.set_expression(y(i) + alpha(0)*y(k-i-1));
    computation y_1("[NN]->{y_1[k,i]: 1<=k<NN and 0<=i<k}", expr(), true, p_float64, global::get_implicit_function());
    y_1.set_expression(z(0,i));
    computation y_2("y_2", {k}, alpha(0));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    y_0.then(beta_0, d)
       .then(alpha_0, d)
       .then(beta_1, computation::root)
       .then(sum_1, k)
       .then(sum_2, k)
       .then(alpha_1, k)
       .then(z, k)
       .then(y_1, k)
       .then(y_2, k);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_r("b_r", {N}, p_float64, a_input);    
    buffer b_y("b_y", {N}, p_float64, a_output);    
    buffer b_z("b_z", {N}, p_float64, a_temporary);    
    buffer b_alpha("b_alpha", {1}, p_float64, a_temporary);    
    buffer b_beta("b_beta", {1}, p_float64, a_temporary);    
    buffer b_sum("b_sum", {1}, p_float64, a_temporary);    

    //Store inputs
    y_0.store_in(&b_y, {});
    beta_0.store_in(&b_beta, {});
    alpha_0.store_in(&b_alpha, {});
    r.store_in(&b_r);  
    y.store_in(&b_y);  
    alpha.store_in(&b_alpha);  
    beta.store_in(&b_beta);
    sum.store_in(&b_sum);

    //Store computations
    beta_1.store_in(&b_beta, {});
    sum_1.store_in(&b_sum, {});
    sum_2.store_in(&b_sum, {});
    alpha_1.store_in(&b_alpha, {});
    z.store_in(&b_z, {i});
    y_1.store_in(&b_y, {i});
    y_2.store_in(&b_y, {k});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_r, &b_y}, "generated_durbin.o");

    return 0;
}
