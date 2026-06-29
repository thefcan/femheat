// M1 (b): assembly of the global system on a tiny 2-element mesh.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "femheat/assembler.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"
#include "femheat/point.hpp"

using femheat::Assembler;
using femheat::Material;
using femheat::Mesh1D;
using femheat::Point;

TEST(Assembler, TwoElementMesh) {
  // Mesh [0, 1] with 2 elements -> h = 0.5, k = 1, so each
  // Ke = (1/0.5) * [[1,-1],[-1,1]] = [[2,-2],[-2,2]].
  const Mesh1D mesh(0.0, 1.0, 2);
  const auto elements = mesh.buildElements();

  // Constant source f = 8 -> per-node element load f*h/2 = 2.
  const auto system = Assembler::assemble(mesh.numNodes(), elements, Material(1.0),
                                          [](const Point&) { return 8.0; });

  // Global stiffness: overlap at the shared middle node sums to 4.
  const Eigen::MatrixXd K(system.K);
  Eigen::MatrixXd expectedK(3, 3);
  expectedK << 2, -2, 0, -2, 4, -2, 0, -2, 2;
  EXPECT_TRUE(K.isApprox(expectedK)) << "K =\n" << K;

  // Global load: middle node receives a contribution from both elements.
  Eigen::VectorXd expectedF(3);
  expectedF << 2, 4, 2;
  EXPECT_TRUE(system.f.isApprox(expectedF)) << "f =\n" << system.f;
}
