#include "femheat/build_info.hpp"

// Including a HighFive header here verifies that the HighFive wrapper compiles
// against the discovered HDF5 installation; calling an HDF5 C function below
// verifies that the HDF5 library actually links.
#include <hdf5.h>

#include <highfive/H5File.hpp>
#include <string>

namespace femheat {

const char* version() {
  return FEMHEAT_VERSION;
}

std::string hdf5Version() {
  unsigned major = 0;
  unsigned minor = 0;
  unsigned release = 0;
  H5get_libversion(&major, &minor, &release);
  return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(release);
}

}  // namespace femheat
