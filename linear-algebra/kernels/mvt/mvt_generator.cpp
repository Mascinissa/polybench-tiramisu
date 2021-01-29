#include <tiramisu/tiramisu.h>
#include "polybench-tiramisu.h"
#include "mvt.h"


using namespace tiramisu;

/*
Matrix vector multiplication composed with another matrix vector multiplication but with transposed matrix.
It takes the following as inputs,
  • A: NxN matrix
  • y1, y2: vectors of length N
and gives the following as outputs:
  • x1: vector of length N, where x1 = x1 + Ay1
  • x2: vector of length N, where x2 = x2 + A^T y2
*/

int main(int argc, char **argv) {
  tiramisu::init("mvt");

  // -------------------------------------------------------
  // Layer I
  // -------------------------------------------------------

  // Iteration variables
  var i("i", 0, N), j("j", 0, N);

  // inputs
  input A("A", {i, j}, p_float64);
  input y1("y1", {j}, p_float64);
  input y2("y2", {j}, p_float64);

  // Computations
  computation x1("x1", {i, j}, p_float64);
  x1.set_expression(x1(i, j) + A(i, j) * y1(j));
  computation x2("x2", {i, j}, p_float64);
  x2.set_expression(x2(i, j) + A(j, i) * y2(j));


  // -------------------------------------------------------
  // Layer II
  // -------------------------------------------------------
  x1.then(x2, computation::root);

  // -------------------------------------------------------
  // Layer III
  // -------------------------------------------------------
  // Input Buffers
  buffer b_A("b_A", {N, N}, p_float64, a_input);
  buffer b_y1("b_y1", {N}, p_float64, a_input);
  buffer b_y2("b_y2", {N}, p_float64, a_input);
  buffer b_x1("b_x1", {N}, p_float64, a_output);
  buffer b_x2("b_x2", {N}, p_float64, a_output);

  // Store inputs
  A.store_in(&b_A);
  y1.store_in(&b_y1);
  y2.store_in(&b_y2);


  // Store computations
  x1.store_in(&b_x1, {i});
  x2.store_in(&b_x2, {i});

  // -------------------------------------------------------
  // Code Generation
  // -------------------------------------------------------
  tiramisu::codegen({&b_A, &b_y1, &b_y2, &b_x1, &b_x2}, "generated_mvt.o");

  return 0;
}
