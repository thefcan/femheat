#include "femheat/mms.hpp"

#include <cmath>

namespace femheat::mms {
namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

double exact1D(double x) { return std::sin(kPi * x); }

double source1D(double x, double k) {
  return k * kPi * kPi * std::sin(kPi * x);
}

double l2Error1D(const Mesh1D& mesh, const Eigen::VectorXd& uh,
                 const std::function<double(double)>& exact) {
  // 3-point Gauss-Legendre nodes/weights on the reference interval [-1, 1].
  const double gp[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
  const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  double sumSq = 0.0;
  for (int e = 0; e < mesh.numElements(); ++e) {
    const double x0 = mesh.nodeX(e);
    const double x1 = mesh.nodeX(e + 1);
    const double h = x1 - x0;
    const double u0 = uh(e);
    const double u1 = uh(e + 1);
    for (int q = 0; q < 3; ++q) {
      const double t = 0.5 * (gp[q] + 1.0);        // map [-1,1] -> [0,1]
      const double x = x0 + t * h;                 // physical point
      const double uhx = u0 * (1.0 - t) + u1 * t;  // linear interpolation
      const double diff = uhx - exact(x);
      sumSq += (h / 2.0) * gw[q] * diff * diff;
    }
  }
  return std::sqrt(sumSq);
}

double observedOrder(double hCoarse, double errCoarse, double hFine,
                     double errFine) {
  return std::log(errCoarse / errFine) / std::log(hCoarse / hFine);
}

}  // namespace femheat::mms
