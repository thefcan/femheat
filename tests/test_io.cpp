// M4: smoke test for the HDF5 + XDMF writer.

#include <gtest/gtest.h>

#include <highfive/H5File.hpp>

#include <Eigen/Dense>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>

#include "femheat/hdf5_writer.hpp"
#include "femheat/tri_mesh.hpp"

using femheat::Hdf5Writer;
using femheat::TriMesh;

TEST(Hdf5Writer, WritesExpectedDatasetsAndXdmf) {
  const TriMesh mesh(1.0, 1.0, 5, 5);  // 25 nodes, 32 triangles
  const Eigen::VectorXd field =
      Eigen::VectorXd::LinSpaced(mesh.numNodes(), 0.0, 1.0);

  const std::string h5 = "femheat_io_test.h5";
  const std::string xdmf = "femheat_io_test.xdmf";
  Hdf5Writer::writeTriMesh(h5, mesh, field);

  // (1) HDF5 datasets exist with the expected shapes.
  {
    HighFive::File file(h5, HighFive::File::ReadOnly);
    ASSERT_TRUE(file.exist("/mesh/coordinates"));
    ASSERT_TRUE(file.exist("/mesh/connectivity"));
    ASSERT_TRUE(file.exist("/solution/T"));

    const auto coords =
        file.getDataSet("/mesh/coordinates").getSpace().getDimensions();
    ASSERT_EQ(coords.size(), 2u);
    EXPECT_EQ(coords[0], static_cast<std::size_t>(mesh.numNodes()));
    EXPECT_EQ(coords[1], 2u);

    const auto conn =
        file.getDataSet("/mesh/connectivity").getSpace().getDimensions();
    ASSERT_EQ(conn.size(), 2u);
    EXPECT_EQ(conn[0], static_cast<std::size_t>(mesh.numElements()));
    EXPECT_EQ(conn[1], 3u);

    const auto t =
        file.getDataSet("/solution/T").getSpace().getDimensions();
    ASSERT_EQ(t.size(), 1u);
    EXPECT_EQ(t[0], static_cast<std::size_t>(mesh.numNodes()));
  }

  // (2) The companion XDMF exists and references the HDF5 datasets.
  {
    std::ifstream in(xdmf);
    ASSERT_TRUE(in.good());
    const std::string content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("TopologyType=\"Triangle\""), std::string::npos);
    EXPECT_NE(content.find("femheat_io_test.h5:/mesh/connectivity"),
              std::string::npos);
    EXPECT_NE(content.find(":/solution/T"), std::string::npos);
  }

  std::remove(h5.c_str());
  std::remove(xdmf.c_str());
}
