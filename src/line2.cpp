#include "femheat/line2.hpp"

#include <cmath>

namespace femheat {

Line2::Line2(int node0, int node1, double x0, double x1)
    : node0_(node0), node1_(node1), x0_(x0), x1_(x1) {}

double Line2::length() const { return std::abs(x1_ - x0_); }

Eigen::MatrixXd Line2::elementStiffness(const Material& material) const {
  const double k = material.conductivity();
  const double c = k / length();  // k / h
  // Ke = (k / h) * [[ 1, -1], [-1, 1]]
  Eigen::MatrixXd ke(2, 2);
  ke << c, -c,
      -c, c;
  return ke;
}

Eigen::VectorXd Line2::elementLoad(double source) const {
  const double v = source * length() / 2.0;  // f * h / 2
  // fe = (f * h / 2) * [1, 1]^T
  Eigen::VectorXd fe(2);
  fe << v, v;
  return fe;
}

std::vector<int> Line2::dofs() const { return {node0_, node1_}; }

Point Line2::centroid() const { return Point{0.5 * (x0_ + x1_), 0.0}; }

}  // namespace femheat
