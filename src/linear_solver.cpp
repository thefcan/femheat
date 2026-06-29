#include "femheat/linear_solver.hpp"

#include <algorithm>

namespace femheat {

LinearSolver::Result LinearSolver::solve(const Eigen::SparseMatrix<double>& K,
                                         const Eigen::VectorXd& b, double tol, int maxIterations) {
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
  cg.setTolerance(tol);
  cg.setMaxIterations(maxIterations > 0 ? static_cast<Eigen::Index>(maxIterations)
                                        : std::max<Eigen::Index>(1000, 20 * K.cols()));
  cg.compute(K);

  Result result;
  result.x = cg.solve(b);
  result.iterations = static_cast<int>(cg.iterations());
  result.error = cg.error();
  result.converged = (cg.info() == Eigen::Success);
  return result;
}

}  // namespace femheat
