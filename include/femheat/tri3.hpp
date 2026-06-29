#ifndef FEMHEAT_TRI3_HPP
#define FEMHEAT_TRI3_HPP

#include <Eigen/Dense>
#include <array>
#include <vector>

#include "femheat/element.hpp"
#include "femheat/material.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Three-node linear triangle element.
///
/// With geometry coefficients
///   b1=y2-y3, b2=y3-y1, b3=y1-y2,    c1=x3-x2, c2=x1-x3, c3=x2-x1,
///   A = |detJ| / 2,
/// the closed-form stiffness and (constant-source) load are
///   Ke(i,j) = (k / (4A)) (b_i b_j + c_i c_j),
///   fe(i)   = f A / 3.
class Tri3 : public IElement {
 public:
  Tri3(std::array<int, 3> nodes, std::array<Point, 3> coords);

  Eigen::MatrixXd elementStiffness(const Material& material) const override;
  Eigen::VectorXd elementLoad(double source) const override;
  std::vector<int> dofs() const override;
  Point centroid() const override;

  double area() const;

 private:
  std::array<int, 3> nodes_;
  std::array<Point, 3> coords_;
};

}  // namespace femheat

#endif  // FEMHEAT_TRI3_HPP
