#ifndef FEMHEAT_MESH1D_HPP
#define FEMHEAT_MESH1D_HPP

#include <memory>
#include <vector>

#include "femheat/element.hpp"

namespace femheat {

/// Uniform 1D mesh of `numElements` two-node segments on [a, b].
/// Nodes are numbered 0..numElements from left to right.
class Mesh1D {
 public:
  Mesh1D(double a, double b, int numElements);

  int numNodes() const { return static_cast<int>(coords_.size()); }
  int numElements() const { return numNodes() - 1; }

  double nodeX(int i) const { return coords_[static_cast<std::size_t>(i)]; }
  const std::vector<double>& coords() const { return coords_; }

  int leftNode() const { return 0; }
  int rightNode() const { return numNodes() - 1; }

  /// Build one Line2 per cell, returned as polymorphic IElement objects.
  std::vector<std::unique_ptr<IElement>> buildElements() const;

 private:
  std::vector<double> coords_;
};

}  // namespace femheat

#endif  // FEMHEAT_MESH1D_HPP
