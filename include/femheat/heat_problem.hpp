#ifndef FEMHEAT_HEAT_PROBLEM_HPP
#define FEMHEAT_HEAT_PROBLEM_HPP

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <vector>

#include "femheat/boundary_condition.hpp"
#include "femheat/element.hpp"
#include "femheat/linear_solver.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"
#include "femheat/point.hpp"
#include "femheat/tri_mesh.hpp"

namespace femheat {

/// Orchestrates a steady-state conduction problem in 1D or 2D. It owns the
/// element set (built from a mesh), the material, the source term and the
/// boundary conditions, then assembles and solves
///   -div( k grad T ) = f
/// returning the nodal temperatures. The same code path serves both the 1D
/// Line2 mesh and the 2D triangle mesh through the IElement interface.
class HeatProblem {
 public:
  HeatProblem(const Mesh1D& mesh, Material material);
  HeatProblem(const TriMesh& mesh, Material material);

  /// Set the volumetric source f as a function of position.
  void setSource(std::function<double(const Point&)> source);

  /// Convenience for a spatially constant source.
  void setConstantSource(double value);

  /// Prescribe T = value at the given node (penalty Dirichlet).
  void addDirichlet(int node, double value);

  /// Add a prescribed boundary flux at the given node (Neumann).
  void addNeumann(int node, double flux);

  /// Assemble, apply boundary conditions and solve. Returns nodal T.
  Eigen::VectorXd solve() const;

  /// Same as solve() but also exposes the solver diagnostics.
  LinearSolver::Result solveDetailed() const;

  int numNodes() const { return numNodes_; }

 private:
  HeatProblem(int numNodes, std::vector<std::unique_ptr<IElement>> elements,
              Material material);

  int numNodes_;
  std::vector<std::unique_ptr<IElement>> elements_;
  Material material_;
  std::function<double(const Point&)> source_;
  std::vector<std::unique_ptr<IBoundaryCondition>> bcs_;
};

}  // namespace femheat

#endif  // FEMHEAT_HEAT_PROBLEM_HPP
