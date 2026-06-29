#ifndef FEMHEAT_LINEAR_SOLVER_HPP
#define FEMHEAT_LINEAR_SOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace femheat {

/// Thin wrapper over Eigen's Conjugate Gradient for the symmetric
/// positive-definite systems produced by the FEM assembly. Eigen's default
/// diagonal (Jacobi) preconditioner keeps the large penalty-BC entries well
/// scaled, so the solve stays accurate.
class LinearSolver {
 public:
  struct Result {
    Eigen::VectorXd x;
    int iterations{0};
    double error{0.0};  // final relative residual reported by Eigen
    bool converged{false};
  };

  /// Solve K x = b. `tol` is the relative residual tolerance; pass
  /// `maxIterations <= 0` to use a generous default.
  static Result solve(const Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& b,
                      double tol = 1.0e-14, int maxIterations = -1);
};

}  // namespace femheat

#endif  // FEMHEAT_LINEAR_SOLVER_HPP
