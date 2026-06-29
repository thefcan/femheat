// M2: 1D Method of Manufactured Solutions — verify second-order convergence.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdio>
#include <vector>

#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"
#include "femheat/mms.hpp"
#include "femheat/point.hpp"

using femheat::HeatProblem;
using femheat::Material;
using femheat::Mesh1D;
using femheat::Point;
namespace mms = femheat::mms;

TEST(MMS1D, SecondOrderConvergence) {
  const double k = 1.0;
  const std::vector<int> ns = {10, 20, 40, 80};
  std::vector<double> hs;
  std::vector<double> errs;

  for (int n : ns) {
    Mesh1D mesh(0.0, 1.0, n);
    HeatProblem problem(mesh, Material(k));
    problem.setSource([k](const Point& p) { return mms::source1D(p.x, k); });
    problem.addDirichlet(mesh.leftNode(), 0.0);
    problem.addDirichlet(mesh.rightNode(), 0.0);

    const Eigen::VectorXd u = problem.solve();
    const double err = mms::l2Error1D(mesh, u, [](double x) { return mms::exact1D(x); });
    hs.push_back(1.0 / n);
    errs.push_back(err);
  }

  // Print the convergence table so the real numbers show up in ctest output.
  std::printf("\n1D MMS convergence (u = sin(pi x)):\n");
  std::printf("  %6s  %12s  %7s\n", "n", "L2 error", "order");
  std::printf("  %6d  %12.4e  %7s\n", ns[0], errs[0], "-");
  for (std::size_t i = 1; i < ns.size(); ++i) {
    const double p = mms::observedOrder(hs[i - 1], errs[i - 1], hs[i], errs[i]);
    std::printf("  %6d  %12.4e  %7.3f\n", ns[i], errs[i], p);
    EXPECT_GE(p, 1.7) << "order n=" << ns[i - 1] << "->" << ns[i] << " = " << p;
    EXPECT_LE(p, 2.3) << "order n=" << ns[i - 1] << "->" << ns[i] << " = " << p;
  }
}
