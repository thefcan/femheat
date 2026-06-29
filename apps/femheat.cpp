// femheat — command-line front-end for the FEM heat-conduction solver.
//
//   femheat --case mms2d --nx 40 --ny 40 --k 1.0 --out result.h5
//
// The mms2d case solves the 2D manufactured solution on the unit square and
// writes result.h5 + result.xdmf for ParaView; mms1d runs the 1D case. With
// --csv it also dumps the node field as x,y,T for quick plotting.

#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include "femheat/build_info.hpp"
#include "femheat/hdf5_writer.hpp"
#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/mesh1d.hpp"
#include "femheat/mms.hpp"
#include "femheat/point.hpp"
#include "femheat/tri_mesh.hpp"

namespace po = boost::program_options;

using femheat::Hdf5Writer;
using femheat::HeatProblem;
using femheat::Material;
using femheat::Mesh1D;
using femheat::Point;
using femheat::TriMesh;
namespace mms = femheat::mms;

namespace {

void writeCsv(const std::string& path, const std::vector<Point>& nodes,
              const Eigen::VectorXd& field) {
  std::ofstream csv(path);
  if (!csv) {
    throw std::runtime_error("femheat: cannot open CSV file: " + path);
  }
  csv << "x,y,T\n";
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    csv << nodes[i].x << "," << nodes[i].y << "," << field(static_cast<Eigen::Index>(i)) << "\n";
  }
}

int runMms2D(int nx, int ny, double k, const std::string& out, const std::string& csv) {
  TriMesh mesh(1.0, 1.0, nx, ny);
  HeatProblem problem(mesh, Material(k));
  problem.setSource([k](const Point& p) { return mms::source2D(p.x, p.y, k); });
  for (int node : mesh.boundaryNodes()) {
    problem.addDirichlet(node, 0.0);
  }

  const Eigen::VectorXd u = problem.solve();
  const double err = mms::l2Error2D(mesh, u, [](double x, double y) { return mms::exact2D(x, y); });
  Hdf5Writer::writeTriMesh(out, mesh, u);
  if (!csv.empty()) {
    writeCsv(csv, mesh.nodes(), u);
  }

  std::cout << "case      : mms2d  (u = sin(pi x) sin(pi y))\n"
            << "mesh      : " << nx << " x " << ny << " nodes, " << mesh.numElements()
            << " triangles\n"
            << "k         : " << k << "\n"
            << "L2 error  : " << err << "\n"
            << "output    : " << out << "  (+ companion .xdmf)\n";
  if (!csv.empty()) {
    std::cout << "csv       : " << csv << "\n";
  }
  return 0;
}

int runMms1D(int nx, double k, const std::string& csv) {
  Mesh1D mesh(0.0, 1.0, nx);
  HeatProblem problem(mesh, Material(k));
  problem.setSource([k](const Point& p) { return mms::source1D(p.x, k); });
  problem.addDirichlet(mesh.leftNode(), 0.0);
  problem.addDirichlet(mesh.rightNode(), 0.0);

  const Eigen::VectorXd u = problem.solve();
  const double err = mms::l2Error1D(mesh, u, [](double x) { return mms::exact1D(x); });

  if (!csv.empty()) {
    std::vector<Point> nodes(static_cast<std::size_t>(mesh.numNodes()));
    for (int i = 0; i < mesh.numNodes(); ++i) {
      nodes[static_cast<std::size_t>(i)] = Point{mesh.nodeX(i), 0.0};
    }
    writeCsv(csv, nodes, u);
  }

  std::cout << "case      : mms1d  (u = sin(pi x))\n"
            << "mesh      : " << nx << " elements\n"
            << "k         : " << k << "\n"
            << "L2 error  : " << err << "\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  po::options_description desc("femheat — 2D finite element heat-conduction solver\nOptions");
  desc.add_options()                              //
      ("help,h", "show this help message")        //
      ("version,v", "print version information")  //
      ("case,c", po::value<std::string>()->default_value("mms2d"),
       "test case: mms1d | mms2d")  //
      ("nx", po::value<int>()->default_value(40),
       "nodes per side (mms2d) / number of elements (mms1d)")  //
      ("ny", po::value<int>()->default_value(40),
       "nodes in y (mms2d only)")  //
      ("k", po::value<double>()->default_value(1.0),
       "thermal conductivity")  //
      ("out,o", po::value<std::string>()->default_value("result.h5"),
       "HDF5 output path (mms2d)")  //
      ("csv", po::value<std::string>()->default_value(""),
       "also write the node field as x,y,T CSV to this path");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "error: " << e.what() << "\n\n" << desc << "\n";
    return 2;
  }

  if (vm.count("help") != 0U) {
    std::cout << desc << "\n";
    return 0;
  }
  if (vm.count("version") != 0U) {
    std::cout << "femheat " << femheat::version() << "\n"
              << "HDF5    " << femheat::hdf5Version() << "\n";
    return 0;
  }

  const auto caseName = vm["case"].as<std::string>();
  const auto nx = vm["nx"].as<int>();
  const auto ny = vm["ny"].as<int>();
  const auto k = vm["k"].as<double>();
  const auto out = vm["out"].as<std::string>();
  const auto csv = vm["csv"].as<std::string>();

  try {
    if (caseName == "mms2d") {
      return runMms2D(nx, ny, k, out, csv);
    }
    if (caseName == "mms1d") {
      return runMms1D(nx, k, csv);
    }
    std::cerr << "error: unknown case '" << caseName << "' (expected mms1d or mms2d)\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
