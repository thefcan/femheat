#include "femheat/tri_mesh.hpp"

#include "femheat/tri3.hpp"

namespace femheat {

TriMesh::TriMesh(double Lx, double Ly, int nx, int ny) : nx_(nx), ny_(ny) {
  const double hx = Lx / (nx - 1);
  const double hy = Ly / (ny - 1);

  nodes_.reserve(static_cast<std::size_t>(nx) * ny);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      nodes_.push_back(Point{i * hx, j * hy});
    }
  }

  tris_.reserve(static_cast<std::size_t>(2) * (nx - 1) * (ny - 1));
  for (int j = 0; j < ny - 1; ++j) {
    for (int i = 0; i < nx - 1; ++i) {
      const int n00 = j * nx_ + i;
      const int n10 = j * nx_ + (i + 1);
      const int n01 = (j + 1) * nx_ + i;
      const int n11 = (j + 1) * nx_ + (i + 1);
      // Two counter-clockwise triangles per cell (positive area).
      tris_.push_back({n00, n10, n11});
      tris_.push_back({n00, n11, n01});
    }
  }
}

std::vector<int> TriMesh::boundaryNodes() const {
  std::vector<int> result;
  for (int j = 0; j < ny_; ++j) {
    for (int i = 0; i < nx_; ++i) {
      if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1) {
        result.push_back(j * nx_ + i);
      }
    }
  }
  return result;
}

std::vector<std::unique_ptr<IElement>> TriMesh::buildElements() const {
  std::vector<std::unique_ptr<IElement>> elements;
  elements.reserve(tris_.size());
  for (const auto& t : tris_) {
    const std::array<Point, 3> coords{nodes_[static_cast<std::size_t>(t[0])],
                                      nodes_[static_cast<std::size_t>(t[1])],
                                      nodes_[static_cast<std::size_t>(t[2])]};
    elements.push_back(std::make_unique<Tri3>(t, coords));
  }
  return elements;
}

}  // namespace femheat
