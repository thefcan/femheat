#include "femheat/heat_problem.hpp"

#include <utility>

#include "femheat/assembler.hpp"

namespace femheat {

HeatProblem::HeatProblem(int numNodes,
                         std::vector<std::unique_ptr<IElement>> elements,
                         Material material)
    : numNodes_(numNodes),
      elements_(std::move(elements)),
      material_(material),
      source_([](const Point&) { return 0.0; }) {}

HeatProblem::HeatProblem(const Mesh1D& mesh, Material material)
    : HeatProblem(mesh.numNodes(), mesh.buildElements(), material) {}

HeatProblem::HeatProblem(const TriMesh& mesh, Material material)
    : HeatProblem(mesh.numNodes(), mesh.buildElements(), material) {}

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
  LinearSystem system =
      Assembler::assemble(numNodes_, elements_, material_, source_);
  for (const auto& bc : bcs_) {
    bc->apply(system.K, system.f);
  }
  return LinearSolver::solve(system.K, system.f);
}

Eigen::VectorXd HeatProblem::solve() const { return solveDetailed().x; }

}  // namespace femheat
