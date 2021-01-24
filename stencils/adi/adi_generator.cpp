#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "adi.h"

using namespace tiramisu;

/*
   TODO description
*/

int main(int argc, char **argv)
{
    double DX, DY, DT;
    double B1, B2;
    double mul1, mul2;
    double a, b, c, d, e, f;

    DX = 1.0/(double)N;
    DY = 1.0/(double)N;
    DT = 1.0/(double)TSTEPS;
    B1 = 2.0;
    B2 = 1.0;
    mul1 = B1 * DT / (DX * DX);
    mul2 = B2 * DT / (DY * DY);

    a = -mul1 /  2.0;
    b = 1.0+mul1;
    c = a;
    d = -mul2 / 2.0;
    e = 1.0+mul2;
    f = d;
    
    tiramisu::init("adi");
    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i("i", 1, N-1), j("j", 1, N-1), t("t", 1, TSTEPS+1);
    var i_f("i_f", 0, N), j_f("j_f", 0, N);

    //inputs
    input u("u", {i, j}, p_float64);

    //Computations
    computation v("v", {i_f,j_f}, 1.0);
    computation q("q", {i_f,j_f}, 1.0);
    computation p("p", {i_f,j_f}, 0.0);

    computation p_col("p_col", {t,i,j}, expr(-c) / (p(i, j-1)*a+b));
    computation q_col("q_col", {t,i,j}, (u(j, i-1)*(-d)+u(j, i)*(1.0+2.0*d) - u(j, i+1)*f-q(i, j-1)*a)/(p(i, j-1)*a+b));
    computation v_col("j_col", {t,i,j}, p(i, j) * v(i, j+1) + q(i, j));

    computation p_row("p_row", {t,i,j}, expr(-f) / (p(i, j-1)*d+e));
    computation q_row("q_row", {t,i,j}, (v(j, i-1)*(-a)+v(j, i)*(1.0+2.0*a) - v(j, i+1)*c-q(i, j-1)*d)/(p(i, j-1)*d+e));
    computation u_row("u_row", {t,i,j}, p(i, j) * u(i, j+1) + q(i, j));

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    v.then(q, j_f)
     .then(p, j_f)
     .then(p_col, computation::root)
     .then(q_col, j)
     .then(v_col, i)
     .then(p_row, t)
     .then(q_row, j)
     .then(u_row, i);

    
    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_u("b_u", {N,N}, p_float64, a_output);    
    buffer b_p("b_p", {N,N}, p_float64, a_temporary);    
    buffer b_q("b_q", {N,N}, p_float64, a_temporary);    
    buffer b_v("b_v", {N,N}, p_float64, a_temporary);    
   

    //Store inputs
    u.store_in(&b_u);

    //Store computations
    v.store_in(&b_v);
    q.store_in(&b_q);
    p.store_in(&b_p);
    p_col.store_in(&b_p,{i,j});
    q_col.store_in(&b_q,{i,j});
    v_col.store_in(&b_v,{i,j});
    p_row.store_in(&b_p,{i,j});
    q_row.store_in(&b_q,{i,j});
    u_row.store_in(&b_u,{i,j});

    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_u}, "generated_adi.o");

    return 0;
}
