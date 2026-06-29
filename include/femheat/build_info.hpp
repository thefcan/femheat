#ifndef FEMHEAT_BUILD_INFO_HPP
#define FEMHEAT_BUILD_INFO_HPP

#include <string>

namespace femheat {

/// Returns the femheat project version string (e.g. "0.1.0").
const char* version();

/// Returns the runtime HDF5 library version (e.g. "1.14.3").
///
/// This calls into HDF5 through the HighFive stack, so a non-empty result
/// confirms the HDF5/HighFive file-I/O toolchain is correctly linked.
std::string hdf5Version();

}  // namespace femheat

#endif  // FEMHEAT_BUILD_INFO_HPP
