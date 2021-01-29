#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "atax.h"


using namespace tiramisu;

/*
Computes A^T times Ax.
It takes the following as inputs,
  • A: MxN matrix
  • x: vector of length N
and gives the following as output:
  • y: vector of length N, where y = A^T (Ax)
*/

int main(int argc, char **argv) {
  tiramisu::init("atax");

  // -------------------------------------------------------
  // Layer I
  // -------------------------------------------------------

  // Iteration variables
  var i("i", 0, M), j("j", 0, N);

  // inputs
  input A("A", {i, j}, p_float64);
  input x("x", {j}, p_float64);

  // Computations
  computation Ax_init("Ax_init", {i}, 0.0);
  computation Ax("Ax", {i, j}, p_float64);
  Ax.set_expression(Ax(i, j) + A(i, j) * x(j));
  computation y_init("y_init", {j}, 0.0);
  computation y("y", {i, j}, p_float64);
  y.set_expression(y(i, j) + A(i, j) * Ax(i, 0));

  // -------------------------------------------------------
  // Layer II
  // -------------------------------------------------------
  Ax_init.then(Ax, computation::root)
         .then(y_init, computation::root)
         .then(y, computation::root);

  // -------------------------------------------------------
  // Layer III
  // -------------------------------------------------------
  // Input Buffers
  buffer b_A("b_A", {M, N}, p_float64, a_input);
  buffer b_Ax("b_Ax", {M}, p_float64, a_temporary);
  buffer b_x("b_x", {N}, p_float64, a_input);
  buffer b_y("b_y", {N}, p_float64, a_output);

  // Store inputs
  A.store_in(&b_A);
  x.store_in(&b_x);

  // Store computations
  Ax_init.store_in(&b_Ax);
  Ax.store_in(&b_Ax, {i});
  y_init.store_in(&b_y);
  y.store_in(&b_y, {j});

  // -------------------------------------------------------
  // Code Generation
  // -------------------------------------------------------
  tiramisu::codegen({&b_A, &b_x, &b_y}, "generated_atax.o");

  return 0;
}
