#ifndef FEMHEAT_HEAT_PROBLEM_HPP
#define FEMHEAT_HEAT_PROBLEM_HPP

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <vector>

#include "femheat/boundary_condition.hpp"
#include "femheat/linear_solver.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"
#include "femheat/point.hpp"

namespace femheat {

/// Orchestrates a 1D steady-state conduction problem. It owns the mesh,
/// material, source term and boundary conditions, then assembles and solves
///   -d/dx( k dT/dx ) = f
/// returning the nodal temperatures.
class HeatProblem {
 public:
  HeatProblem(Mesh1D mesh, Material material);

  /// Set the volumetric source f as a function of position.
  void setSource(std::function<double(const Point&)> source);

  /// Convenience for a spatially constant source.
  void setConstantSource(double value);

  /// Prescribe T = value at the given node (penalty Dirichlet).
  void addDirichlet(int node, double value);

  /// Assemble, apply boundary conditions and solve. Returns nodal T.
  Eigen::VectorXd solve() const;

  /// Same as solve() but also exposes the solver diagnostics.
  LinearSolver::Result solveDetailed() const;

  const Mesh1D& mesh() const { return mesh_; }

 private:
  Mesh1D mesh_;
  Material material_;
  std::function<double(const Point&)> source_;
  std::vector<std::unique_ptr<IBoundaryCondition>> bcs_;
};

}  // namespace femheat

#endif  // FEMHEAT_HEAT_PROBLEM_HPP
