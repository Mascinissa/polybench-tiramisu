#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "bicg.h"


using namespace tiramisu;

/*
Kernel of BiCGSTAB (BiConjugate Gradient STABilized method).
It takes the following as inputs,
  • A: NxM matrix
  • p: vector of length M
  • r: vector of length N
and gives the following as output:
  • q: vector of length N, where q = A*p
  • s: vector of length M, where s = A^T*r
*/

int main(int argc, char **argv) {
  tiramisu::init("bicg");

  // -------------------------------------------------------
  // Layer I
  // -------------------------------------------------------

  // Iteration variables
  var i("i", 0, N), j("j", 0, M);

  // inputs
  input A("A", {i, j}, p_float64);
  input p("p", {i}, p_float64);
  input r("r", {j}, p_float64);

  // Computations
  computation q_init("q_init", {i}, 0.0);
  computation q("q", {i, j}, p_float64);
  q.set_expression(q(i, j) + A(i, j) * p(j));
  computation s_init("s_init", {j}, 0.0);
  computation s("s", {i, j}, p_float64);
  s.set_expression(s(i, j) + A(i, j) * r(i));

  // -------------------------------------------------------
  // Layer II
  // -------------------------------------------------------
  s_init.then(q_init, computation::root)
        .then(s, i)
        .then(q, j);

  // -------------------------------------------------------
  // Layer III
  // -------------------------------------------------------
  // Input Buffers
  buffer b_A("b_A", {N, M}, p_float64, a_input);
  buffer b_p("b_p", {M}, p_float64, a_input);
  buffer b_r("b_r", {N}, p_float64, a_input);
  buffer b_q("b_q", {N}, p_float64, a_output);
  buffer b_s("b_s", {M}, p_float64, a_output);

  // Store inputs
  A.store_in(&b_A);
  p.store_in(&b_p);
  r.store_in(&b_r);
  q.store_in(&b_q);
  s.store_in(&b_s);

  // Store computations
  q_init.store_in(&b_q);
  q.store_in(&b_q, {i});
  s_init.store_in(&b_s);
  s.store_in(&b_s, {j});

  // -------------------------------------------------------
  // Code Generation
  // -------------------------------------------------------
  tiramisu::codegen({&b_A, &b_p, &b_r, &b_q, &b_s}, "generated_bicg.o");

  return 0;
}
