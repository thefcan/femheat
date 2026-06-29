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

double exact2D(double x, double y) {
  return std::sin(kPi * x) * std::sin(kPi * y);
}

double source2D(double x, double y, double k) {
  return 2.0 * k * kPi * kPi * std::sin(kPi * x) * std::sin(kPi * y);
}

double l2Error2D(const TriMesh& mesh, const Eigen::VectorXd& uh,
                 const std::function<double(double, double)>& exact) {
  // 3-point barycentric quadrature (weight 1/3 at each point).
  const double bary[3][3] = {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
                             {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
                             {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}};
  const double w = 1.0 / 3.0;
  const auto& nodes = mesh.nodes();

  double sumSq = 0.0;
  for (const auto& t : mesh.triangles()) {
    const Point& a = nodes[static_cast<std::size_t>(t[0])];
    const Point& b = nodes[static_cast<std::size_t>(t[1])];
    const Point& c = nodes[static_cast<std::size_t>(t[2])];
    const double ua = uh(t[0]);
    const double ub = uh(t[1]);
    const double uc = uh(t[2]);
    const double detJ = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
    const double area = std::abs(detJ) / 2.0;
    for (int q = 0; q < 3; ++q) {
      const double l1 = bary[q][0];
      const double l2 = bary[q][1];
      const double l3 = bary[q][2];
      const double x = l1 * a.x + l2 * b.x + l3 * c.x;
      const double y = l1 * a.y + l2 * b.y + l3 * c.y;
      const double uhq = l1 * ua + l2 * ub + l3 * uc;
      const double diff = uhq - exact(x, y);
      sumSq += area * w * diff * diff;
    }
  }
  return std::sqrt(sumSq);
}

double observedOrder(double hCoarse, double errCoarse, double hFine,
                     double errFine) {
  return std::log(errCoarse / errFine) / std::log(hCoarse / hFine);
}

}  // namespace femheat::mms
