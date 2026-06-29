// M3: element-level unit tests for the 2D linear triangle Tri3.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <array>

#include "femheat/material.hpp"
#include "femheat/point.hpp"
#include "femheat/tri3.hpp"

using femheat::Material;
using femheat::Point;
using femheat::Tri3;

TEST(Tri3, StiffnessReferenceTriangle) {
  // Right triangle (0,0), (1,0), (0,1) with k = 1. The closed-form stiffness
  // is the classic linear-triangle Laplacian.
  const Tri3 e({0, 1, 2}, {Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.0, 1.0}});
  const Eigen::MatrixXd ke = e.elementStiffness(Material(1.0));

  Eigen::MatrixXd expected(3, 3);
  expected << 1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(ke(i, j), expected(i, j), 1e-12) << "ke(" << i << "," << j << ")";
    }
  }

  // A constant temperature produces no flux, so each row sums to zero.
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(ke.row(i).sum(), 0.0, 1e-12);
  }
}

TEST(Tri3, AreaLoadAndCentroid) {
  const Tri3 e({0, 1, 2}, {Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.0, 1.0}});
  EXPECT_DOUBLE_EQ(e.area(), 0.5);

  const Eigen::VectorXd fe = e.elementLoad(6.0);  // f*A/3 = 1 at each node
  EXPECT_DOUBLE_EQ(fe(0), 1.0);
  EXPECT_DOUBLE_EQ(fe(1), 1.0);
  EXPECT_DOUBLE_EQ(fe(2), 1.0);

  EXPECT_DOUBLE_EQ(e.centroid().x, 1.0 / 3.0);
  EXPECT_DOUBLE_EQ(e.centroid().y, 1.0 / 3.0);
}
