#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "doitgen.h"

using namespace tiramisu;

/*
It takes the following as inputs,
  •A:R×Q×S array
  •x:P×S array
and gives the following as output:
  •Aout:R×Q×P array A(r,q,p) =∑_s A(r,q,s)x(p,s)
*/

int main(int argc, char **argv) {
  tiramisu::init("doitgen");

  // -------------------------------------------------------
  // Layer I
  // -------------------------------------------------------

  // Iteration variables
  var r("r", 0, R), q("q", 0, Q), p("p", 0, P), s("s", 0, P);

  // inputs
  input A("A", {r, q, s}, p_float64);
  input x("x", {p, s}, p_float64);
  input sum("sum", {p}, p_float64);

  // Computations
  computation sum_init("sum_init", {r, q, p}, 0.0);
  computation sum_comp("sum_comp", {r, q, p, s}, sum(p) + A(r, q, s) * x(s, p));
  computation A_out("A_out", {r, q, p}, sum(p));


  // -------------------------------------------------------
  // Layer II
  // -------------------------------------------------------
    sum_init.then(sum_comp, p)
            .then(A_out, q);

  // -------------------------------------------------------
  // Layer III
  // -------------------------------------------------------
  // Input Buffers
  buffer b_A("b_A", {R, Q, P}, p_float64, a_output);
  buffer b_sum("b_sum", {P}, p_float64, a_temporary);
  buffer b_x("b_x", {P, P}, p_float64, a_input);

  // Store inputs
  sum.store_in(&b_sum);
  A.store_in(&b_A);
  x.store_in(&b_x);

  // Store computations
  sum_init.store_in(&b_sum, {p});
  sum_comp.store_in(&b_sum, {p});
  A_out.store_in(&b_A, {r, q, p});

  // -------------------------------------------------------
  // Code Generation
  // -------------------------------------------------------
  tiramisu::codegen({&b_A, &b_x}, "generated_doitgen.o");

  return 0;
}
