#include "femheat/boundary_condition.hpp"

namespace femheat {

DirichletBC::DirichletBC(int node, double value, double beta)
    : node_(node), value_(value), beta_(beta) {}

void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& rhs) const {
  // Penalty method: inflate the diagonal and the rhs so node_'s equation
  // collapses to (K_ii + beta) * u_i ≈ beta * value, i.e. u_i ≈ value.
  K.coeffRef(node_, node_) += beta_;
  rhs(node_) += beta_ * value_;
}

NeumannBC::NeumannBC(int node, double flux) : node_(node), flux_(flux) {}

void NeumannBC::apply(Eigen::SparseMatrix<double>& /*K*/, Eigen::VectorXd& rhs) const {
  // A prescribed flux enters the weak form as a natural boundary term, added
  // straight to the load vector — the stiffness matrix is untouched.
  rhs(node_) += flux_;
}

}  // namespace femheat
