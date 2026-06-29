#ifndef FEMHEAT_HDF5_WRITER_HPP
#define FEMHEAT_HDF5_WRITER_HPP

#include <Eigen/Dense>
#include <string>

#include "femheat/tri_mesh.hpp"

namespace femheat {

/// Writes a 2D triangle mesh and a nodal scalar field to an HDF5 file and emits
/// a companion XDMF file (same basename, ".xdmf") so the result opens directly
/// in ParaView.
///
/// HDF5 layout:
///   /mesh/coordinates   double [numNodes][2]
///   /mesh/connectivity  int    [numElements][3]   (0-based node indices)
///   /solution/<field>   double [numNodes]
class Hdf5Writer {
 public:
  static void writeTriMesh(const std::string& h5Path, const TriMesh& mesh,
                           const Eigen::VectorXd& field,
                           const std::string& fieldName = "T");
};

}  // namespace femheat

#endif  // FEMHEAT_HDF5_WRITER_HPP
