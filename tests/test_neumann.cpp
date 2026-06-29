// M2: Neumann (flux) boundary condition on a mixed boundary-value problem.

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"

using femheat::HeatProblem;
using femheat::Material;
using femheat::Mesh1D;

TEST(Neumann, MixedDirichletNeumannMatchesAnalytic) {
  // -k u'' = f on [0, 1], u(0) = 0, k u'(1) = q.
  // Exact: u(x) = -(f / (2k)) x^2 + ((q + f) / k) x.
  const double k = 1.0;
  const double f = 2.0;
  const double q = 0.5;
  const Mesh1D mesh(0.0, 1.0, 10);

  HeatProblem problem(mesh, Material(k));
  problem.setConstantSource(f);
  problem.addDirichlet(mesh.leftNode(), 0.0);
  problem.addNeumann(mesh.rightNode(), q);  // +q at the right end (normal +x)

  const Eigen::VectorXd u = problem.solve();
  for (int i = 0; i < mesh.numNodes(); ++i) {
    const double x = mesh.nodeX(i);
    const double exact = -(f / (2.0 * k)) * x * x + ((q + f) / k) * x;
    EXPECT_NEAR(u(i), exact, 1e-10) << "x = " << x;
  }
}
