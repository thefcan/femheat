#ifndef FEMHEAT_LINE2_HPP
#define FEMHEAT_LINE2_HPP

#include <Eigen/Dense>
#include <vector>

#include "femheat/element.hpp"
#include "femheat/material.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Two-node linear element on a 1D segment.
///
/// For an element of length h and conductivity k:
///   Ke = (k / h) * [[ 1, -1],
///                   [-1,  1]]
///   fe = (f * h / 2) * [1, 1]^T          (constant source f)
class Line2 : public IElement {
 public:
  Line2(int node0, int node1, double x0, double x1);

  Eigen::MatrixXd elementStiffness(const Material& material) const override;
  Eigen::VectorXd elementLoad(double source) const override;
  std::vector<int> dofs() const override;
  Point centroid() const override;

  double length() const;

 private:
  int node0_;
  int node1_;
  double x0_;
  double x1_;
};

}  // namespace femheat

#endif  // FEMHEAT_LINE2_HPP
