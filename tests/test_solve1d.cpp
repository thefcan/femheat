// M1 (c): a full 1D solve must reproduce the analytic solution.
//
// For -k u'' = f (constant f) on [0, L] with u(0)=u0 and u(L)=uL, the exact
// solution is
//   u(x) = u0 + (uL - u0) x / L + (f / (2k)) x (L - x).
// 1D linear FEM is *nodally exact* for this problem, so the computed nodal
// temperatures must match the analytic values to ~machine precision.

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"

using femheat::HeatProblem;
using femheat::Material;
using femheat::Mesh1D;

TEST(Solve1D, ConstantSourceZeroBoundary) {
  const double k = 1.0;
  const double f = 2.0;  // exact: u(x) = x (1 - x)
  const Mesh1D mesh(0.0, 1.0, 10);

  HeatProblem problem(mesh, Material(k));
  problem.setConstantSource(f);
  problem.addDirichlet(mesh.leftNode(), 0.0);
  problem.addDirichlet(mesh.rightNode(), 0.0);

  const auto result = problem.solveDetailed();
  EXPECT_TRUE(result.converged);

  for (int i = 0; i < mesh.numNodes(); ++i) {
    const double x = mesh.nodeX(i);
    const double exact = x * (1.0 - x);
    EXPECT_NEAR(result.x(i), exact, 1e-10) << "node " << i << " x = " << x;
  }
}

TEST(Solve1D, ConstantSourceGeneralBoundary) {
  const double k = 2.0;
  const double f = 4.0;
  const double u0 = 1.0;
  const double uL = 3.0;
  const double L = 1.0;
  const Mesh1D mesh(0.0, L, 5);  // k/h = 10 keeps the penalty error ~1e-11

  HeatProblem problem(mesh, Material(k));
  problem.setConstantSource(f);
  problem.addDirichlet(mesh.leftNode(), u0);
  problem.addDirichlet(mesh.rightNode(), uL);

  const Eigen::VectorXd u = problem.solve();
  for (int i = 0; i < mesh.numNodes(); ++i) {
    const double x = mesh.nodeX(i);
    const double exact = u0 + (uL - u0) * x / L + (f / (2.0 * k)) * x * (L - x);
    EXPECT_NEAR(u(i), exact, 1e-10) << "node " << i << " x = " << x;
  }
}
