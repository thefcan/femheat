#include "femheat/heat_problem.hpp"

#include <utility>

#include "femheat/assembler.hpp"

namespace femheat {

HeatProblem::HeatProblem(Mesh1D mesh, Material material)
    : mesh_(std::move(mesh)),
      material_(material),
      source_([](const Point&) { return 0.0; }) {}

void HeatProblem::setSource(std::function<double(const Point&)> source) {
  source_ = std::move(source);
}

void HeatProblem::setConstantSource(double value) {
  source_ = [value](const Point&) { return value; };
}

void HeatProblem::addDirichlet(int node, double value) {
  bcs_.push_back(std::make_unique<DirichletBC>(node, value));
}

void HeatProblem::addNeumann(int node, double flux) {
  bcs_.push_back(std::make_unique<NeumannBC>(node, flux));
}

LinearSolver::Result HeatProblem::solveDetailed() const {
  const auto elements = mesh_.buildElements();
  LinearSystem system =
      Assembler::assemble(mesh_.numNodes(), elements, material_, source_);
  for (const auto& bc : bcs_) {
    bc->apply(system.K, system.f);
  }
  return LinearSolver::solve(system.K, system.f);
}

Eigen::VectorXd HeatProblem::solve() const { return solveDetailed().x; }

}  // namespace femheat
