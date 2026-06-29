#ifndef FEMHEAT_BOUNDARY_CONDITION_HPP
#define FEMHEAT_BOUNDARY_CONDITION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace femheat {

/// Strategy interface for boundary conditions: each implementation knows how to
/// modify the assembled system (K, rhs) so that it enforces itself.
class IBoundaryCondition {
 public:
  virtual ~IBoundaryCondition() = default;
  virtual void apply(Eigen::SparseMatrix<double>& K,
                     Eigen::VectorXd& rhs) const = 0;
};

/// Dirichlet condition T = value at a single node, enforced with the penalty
/// method:
///   K(i,i) += beta;   rhs(i) += beta * value.
/// A large beta makes node i's equation dominate, so the solution there tends
/// to `value` while keeping the matrix symmetric, sparse and positive definite.
/// Eigen's diagonal (Jacobi) preconditioner rescales the inflated diagonal, so
/// the conjugate-gradient solve stays accurate despite the large penalty.
class DirichletBC : public IBoundaryCondition {
 public:
  static constexpr double kDefaultPenalty = 1.0e12;

  DirichletBC(int node, double value, double beta = kDefaultPenalty);

  void apply(Eigen::SparseMatrix<double>& K,
             Eigen::VectorXd& rhs) const override;

 private:
  int node_;
  double value_;
  double beta_;
};

}  // namespace femheat

#endif  // FEMHEAT_BOUNDARY_CONDITION_HPP
