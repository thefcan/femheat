// M3: 2D Method of Manufactured Solutions — verify second-order convergence
// of the linear-triangle solver on the unit square.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdio>
#include <vector>

#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/mms.hpp"
#include "femheat/point.hpp"
#include "femheat/tri_mesh.hpp"

using femheat::HeatProblem;
using femheat::Material;
using femheat::Point;
using femheat::TriMesh;
namespace mms = femheat::mms;

TEST(MMS2D, SecondOrderConvergence) {
  const double k = 1.0;
  const std::vector<int> ns = {11, 21, 41};  // nodes per side; h = 1/(n-1)
  std::vector<double> hs;
  std::vector<double> errs;

  for (int n : ns) {
    TriMesh mesh(1.0, 1.0, n, n);
    HeatProblem problem(mesh, Material(k));
    problem.setSource([k](const Point& p) { return mms::source2D(p.x, p.y, k); });
    for (int node : mesh.boundaryNodes()) {
      problem.addDirichlet(node, 0.0);
    }

    const auto result = problem.solveDetailed();
    EXPECT_LT(result.error, 1e-9) << "CG residual too large for n=" << n;

    const double err = mms::l2Error2D(
        mesh, result.x, [](double x, double y) { return mms::exact2D(x, y); });
    hs.push_back(1.0 / (n - 1));
    errs.push_back(err);
  }

  std::printf("\n2D MMS convergence (u = sin(pi x) sin(pi y)):\n");
  std::printf("  %10s  %12s  %7s\n", "nodes/side", "L2 error", "order");
  std::printf("  %10d  %12.4e  %7s\n", ns[0], errs[0], "-");
  for (std::size_t i = 1; i < ns.size(); ++i) {
    const double p = mms::observedOrder(hs[i - 1], errs[i - 1], hs[i], errs[i]);
    std::printf("  %10d  %12.4e  %7.3f\n", ns[i], errs[i], p);
    EXPECT_GE(p, 1.7) << "order " << ns[i - 1] << "->" << ns[i] << " = " << p;
    EXPECT_LE(p, 2.3) << "order " << ns[i - 1] << "->" << ns[i] << " = " << p;
  }
}
