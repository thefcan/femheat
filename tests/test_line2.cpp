// M1 (a): element-level unit tests for the 1D linear element Line2.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <vector>

#include "femheat/line2.hpp"
#include "femheat/material.hpp"

using femheat::Line2;
using femheat::Material;

TEST(Line2, StiffnessMatrixEntries) {
  // Element on [0, 0.5] (h = 0.5) with k = 2  ->  k/h = 4.
  const Line2 e(0, 1, 0.0, 0.5);
  const Eigen::MatrixXd ke = e.elementStiffness(Material(2.0));

  const double c = 4.0;  // k / h
  EXPECT_DOUBLE_EQ(ke(0, 0), c);
  EXPECT_DOUBLE_EQ(ke(0, 1), -c);
  EXPECT_DOUBLE_EQ(ke(1, 0), -c);
  EXPECT_DOUBLE_EQ(ke(1, 1), c);
}

TEST(Line2, LoadVectorEntries) {
  // h = 0.5, constant source f = 8  ->  f*h/2 = 2 on each node.
  const Line2 e(0, 1, 0.0, 0.5);
  const Eigen::VectorXd fe = e.elementLoad(8.0);
  EXPECT_DOUBLE_EQ(fe(0), 2.0);
  EXPECT_DOUBLE_EQ(fe(1), 2.0);
}

TEST(Line2, DofsLengthAndCentroid) {
  const Line2 e(3, 4, 1.0, 2.0);
  EXPECT_EQ(e.dofs(), (std::vector<int>{3, 4}));
  EXPECT_DOUBLE_EQ(e.length(), 1.0);
  EXPECT_DOUBLE_EQ(e.centroid().x, 1.5);
}
