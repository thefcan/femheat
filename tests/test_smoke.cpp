// M0 smoke test.
//
// Confirms that femheat_lib links and that the third-party toolchain is wired
// up correctly: GoogleTest runs this binary, and the library pulls in
// HDF5/HighFive. The real numerical unit/regression tests arrive in M1+.

#include <gtest/gtest.h>

#include <string>

#include "femheat/build_info.hpp"

TEST(Smoke, ProjectVersionIsReported) {
  const std::string v = femheat::version();
  EXPECT_FALSE(v.empty());
  EXPECT_EQ(v, "0.1.0");
}

TEST(Smoke, Hdf5ToolchainLinks) {
  // A non-empty, non-zero version means HDF5 is reachable through HighFive,
  // i.e. the file-I/O stack used from M4 onward is correctly linked.
  const std::string h = femheat::hdf5Version();
  EXPECT_FALSE(h.empty());
  EXPECT_NE(h, "0.0.0");
}
