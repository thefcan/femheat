#ifndef FEMHEAT_TRI_MESH_HPP
#define FEMHEAT_TRI_MESH_HPP

#include <array>
#include <memory>
#include <vector>

#include "femheat/element.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Structured triangulation of the rectangle [0, Lx] x [0, Ly] using an
/// nx-by-ny node grid. Each grid cell is split into two triangles, giving
/// 2*(nx-1)*(ny-1) elements. Node (i, j) has index j*nx + i.
class TriMesh {
 public:
  TriMesh(double Lx, double Ly, int nx, int ny);

  int numNodes() const { return static_cast<int>(nodes_.size()); }
  int numElements() const { return static_cast<int>(tris_.size()); }

  const std::vector<Point>& nodes() const { return nodes_; }
  const std::vector<std::array<int, 3>>& triangles() const { return tris_; }

  /// Indices of all nodes lying on the rectangle boundary.
  std::vector<int> boundaryNodes() const;

  /// Build the Tri3 elements as polymorphic IElement objects.
  std::vector<std::unique_ptr<IElement>> buildElements() const;

 private:
  int nx_;
  int ny_;
  std::vector<Point> nodes_;
  std::vector<std::array<int, 3>> tris_;
};

}  // namespace femheat

#endif  // FEMHEAT_TRI_MESH_HPP
