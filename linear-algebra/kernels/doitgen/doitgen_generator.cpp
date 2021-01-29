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

  // Computations
  computation A_init("A_init", {r, q, p}, 0.0);
  computation A_sum("A_sum", {r, q, p, s}, p_float64);
  A_sum.set_expression(A_sum(r, q, p, s) + A(r, q, s) * x(p, s));
  computation A_out("A_out", {r, q, p}, A_sum(r, q, p, 0));


  // -------------------------------------------------------
  // Layer II
  // -------------------------------------------------------
  A_init.then(A_sum, computation::root)
        .then(A_out, computation::root);
  // -------------------------------------------------------
  // Layer III
  // -------------------------------------------------------
  // Input Buffers
  buffer b_A("b_A", {R, Q, P}, p_float64, a_output);
  buffer b_A_sum("b_A_sum", {R, Q, P}, p_float64, a_temporary);
  buffer b_x("b_x", {P, P}, p_float64, a_input);


  // Store inputs
  A.store_in(&b_A);
  x.store_in(&b_x);

  // Store computations
  A_init.store_in(&b_A_sum);
  A_sum.store_in(&b_A_sum, {r, q, p});
  A_out.store_in(&b_A, {r, q, p});

  // -------------------------------------------------------
  // Code Generation
  // -------------------------------------------------------
  tiramisu::codegen({&b_A, &b_x}, "generated_doitgen.o");

  return 0;
}
