// femheat — command-line front-end for the FEM heat-conduction solver.
//
// At milestone M0 this only wires up Boost.Program_options and reports build
// information; the actual solver cases (--case mms1d/mms2d ...) are added in
// later milestones.

#include <boost/program_options.hpp>

#include <iostream>

#include "femheat/build_info.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv) {
  po::options_description desc(
      "femheat — 2D finite element heat-conduction solver\nOptions");
  desc.add_options()                               //
      ("help,h", "show this help message")         //
      ("version,v", "print version information");  //

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

  std::cout << "femheat " << femheat::version()
            << " — run with --help for usage.\n"
            << "Solver cases are added in milestones M1-M4.\n";
  return 0;
}
