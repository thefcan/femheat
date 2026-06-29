#ifndef FEMHEAT_ELEMENT_HPP
#define FEMHEAT_ELEMENT_HPP

#include <Eigen/Dense>
#include <vector>

#include "femheat/material.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Abstract finite element. Concrete elements (Line2, Tri3) implement the local
/// stiffness matrix, the local load vector, and report the global node indices
/// (degrees of freedom) they couple. Working through this interface keeps the
/// Assembler independent of the element type — this is the polymorphism that
/// lets the same assembly code serve both the 1D and 2D solvers.
class IElement {
 public:
  virtual ~IElement() = default;

  /// Local stiffness matrix Ke (n_nodes x n_nodes) for the given material.
  virtual Eigen::MatrixXd elementStiffness(const Material& material) const = 0;

  /// Local load vector fe (n_nodes) for a source assumed constant over the
  /// element with value `source`.
  virtual Eigen::VectorXd elementLoad(double source) const = 0;

  /// Global node indices coupled by this element, in local order.
  virtual std::vector<int> dofs() const = 0;

  /// Element centroid, used to sample the source term.
  virtual Point centroid() const = 0;
};

}  // namespace femheat

#endif  // FEMHEAT_ELEMENT_HPP
