#ifndef FEMHEAT_MMS_HPP
#define FEMHEAT_MMS_HPP

#include <Eigen/Dense>
#include <functional>

#include "femheat/mesh1d.hpp"

namespace femheat::mms {

/// 1D manufactured solution on [0, 1]:
///   u(x) = sin(pi x)   ->   f(x) = k pi^2 sin(pi x),   u(0) = u(1) = 0.
double exact1D(double x);
double source1D(double x, double k);

/// Discrete L2 norm of (u_h - exact) over a 1D mesh, where u_h is the FEM
/// nodal solution interpolated linearly inside each element. Each element is
/// integrated with a 3-point Gauss-Legendre rule.
double l2Error1D(const Mesh1D& mesh, const Eigen::VectorXd& uh,
                 const std::function<double(double)>& exact);

/// Observed convergence order from two (mesh size, error) samples:
///   p = log(errCoarse / errFine) / log(hCoarse / hFine).
double observedOrder(double hCoarse, double errCoarse, double hFine,
                     double errFine);

}  // namespace femheat::mms

#endif  // FEMHEAT_MMS_HPP
