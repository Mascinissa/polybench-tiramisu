#include <tiramisu/tiramisu.h>
#include "fdtd_2d.h"
#include "polybench-tiramisu.h"

using namespace tiramisu;

/*
Simplified Finite-Difference Time-Domain method for 2D data models electric and magnetic fields
based on Maxwell’s equations. In particular, the polarization used here is T Ez; Transverse Electric in z
direction. It is a stencil involving three variables, Ex, Ey, and Hz. Ex and Ey are electric fields varying in x
and y axes, where Hz is the magnetic field along z axis.
*/

int main(int argc, char **argv)
{
    tiramisu::init("fdtd_2d");

    // -------------------------------------------------------
    // Layer I
    // ------------------------------------------------------- 

    //Iteration variables    
    var i_f("i_f", 0, NX), j_f("j_f", 0, NY), i_m("i_m", 0, NX-1), j_m("j_m", 0, NY-1);
    var t("t", 0, TMAX), i("i", 1, NX), j("j", 1, NY);
    
    //inputs
    input fict("fict", {t}, p_float64);
    input ey("ey", {i_f, j_f}, p_float64);
    input ex("ex", {i_f, j_f}, p_float64);
    input hz("hz", {i_f, j_f}, p_float64);

    //Computations
    computation ey_slice("ey_slice", {t,j_f}, fict(t));
    computation ey_out("ey_out", {t, i, j_f}, ey(i, j_f) - (hz(i, j_f) - hz(i-1, j_f))*0.5);
    computation ex_out("ex_out", {t, i_f, j}, ex(i_f, j) - (hz(i_f, j) - hz(i_f, j - 1))*0.5);
    computation hz_out("hz_out", {t, i_m, j_m}, hz(i_m, j_m) - (ex(i_m, j_m + 1) - ex(i_m, j_m) + ey(i_m + 1, j_m) - ey(i_m, j_m))*0.7);

    // -------------------------------------------------------
    // Layer II
    // -------------------------------------------------------
    ey_slice.then(ey_out, t)
            .then(ex_out, t)
            .then(hz_out, t);

    // -------------------------------------------------------
    // Layer III
    // -------------------------------------------------------
    //Input Buffers
    buffer b_fict("b_fict", {TMAX}, p_float64, a_input);    
    buffer b_ey("b_ey", {NX,NY}, p_float64, a_output);    
    buffer b_ex("b_ex", {NX,NY}, p_float64, a_output);    
    buffer b_hz("b_hz", {NX,NY}, p_float64, a_output);    

    //Store inputs
    fict.store_in(&b_fict);
    ey.store_in(&b_ey);
    ex.store_in(&b_ex);
    hz.store_in(&b_hz);


    //Store computations
    ey_slice.set_access("{ey_slice[t,j_f]->b_ey[0,j_f]}");
    ey_out.store_in(&b_ey, {i,j_f});
    ex_out.store_in(&b_ex, {i_f,j});
    hz_out.store_in(&b_hz, {i_m, j_m});


    // -------------------------------------------------------
    // Code Generation
    // -------------------------------------------------------
    tiramisu::codegen({&b_ex, &b_ey, &b_hz, &b_fict}, "generated_fdtd_2d.o");

    return 0;
}
