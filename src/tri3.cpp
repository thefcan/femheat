#include "femheat/tri3.hpp"

#include <cmath>

namespace femheat {

namespace {
// 2 * signed area of the triangle (the Jacobian determinant).
double signedDetJ(const std::array<Point, 3>& p) {
  return (p[1].x - p[0].x) * (p[2].y - p[0].y) -
         (p[2].x - p[0].x) * (p[1].y - p[0].y);
}
}  // namespace

Tri3::Tri3(std::array<int, 3> nodes, std::array<Point, 3> coords)
    : nodes_(nodes), coords_(coords) {}

double Tri3::area() const { return std::abs(signedDetJ(coords_)) / 2.0; }

Eigen::MatrixXd Tri3::elementStiffness(const Material& material) const {
  const Point& p1 = coords_[0];
  const Point& p2 = coords_[1];
  const Point& p3 = coords_[2];

  // Gradients of the linear shape functions are constant over the element.
  const double b[3] = {p2.y - p3.y, p3.y - p1.y, p1.y - p2.y};
  const double c[3] = {p3.x - p2.x, p1.x - p3.x, p2.x - p1.x};

  const double coeff = material.conductivity() / (4.0 * area());
  Eigen::MatrixXd ke(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ke(i, j) = coeff * (b[i] * b[j] + c[i] * c[j]);
    }
  }
  return ke;
}

Eigen::VectorXd Tri3::elementLoad(double source) const {
  // fe(i) = f * A / 3 — a constant source is split equally among the 3 nodes.
  Eigen::VectorXd fe(3);
  fe.setConstant(source * area() / 3.0);
  return fe;
}

std::vector<int> Tri3::dofs() const {
  return {nodes_[0], nodes_[1], nodes_[2]};
}

Point Tri3::centroid() const {
  return Point{(coords_[0].x + coords_[1].x + coords_[2].x) / 3.0,
               (coords_[0].y + coords_[1].y + coords_[2].y) / 3.0};
}

}  // namespace femheat
