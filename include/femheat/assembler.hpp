#ifndef FEMHEAT_ASSEMBLER_HPP
#define FEMHEAT_ASSEMBLER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <memory>
#include <vector>

#include "femheat/element.hpp"
#include "femheat/material.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Global linear system K u = f.
struct LinearSystem {
  Eigen::SparseMatrix<double> K;
  Eigen::VectorXd f;
};

/// Scatters element stiffness matrices and load vectors into a global sparse
/// system. The source term is sampled at each element's centroid (exact for a
/// constant source; second-order accurate for smooth sources).
class Assembler {
 public:
  static LinearSystem assemble(int numNodes, const std::vector<std::unique_ptr<IElement>>& elements,
                               const Material& material,
                               const std::function<double(const Point&)>& source);
};

}  // namespace femheat

#endif  // FEMHEAT_ASSEMBLER_HPP
