#include "femheat/hdf5_writer.hpp"

#include <highfive/H5File.hpp>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace femheat {
namespace {

// File name without any leading directory, used inside the XDMF references.
std::string h5BaseName(const std::string& path) {
  const auto slash = path.find_last_of("/\\");
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

// Companion XDMF path: replace a trailing ".h5" with ".xdmf".
std::string xdmfPathFor(const std::string& h5Path) {
  const std::string suffix = ".h5";
  if (h5Path.size() >= suffix.size() &&
      h5Path.compare(h5Path.size() - suffix.size(), suffix.size(), suffix) ==
          0) {
    return h5Path.substr(0, h5Path.size() - suffix.size()) + ".xdmf";
  }
  return h5Path + ".xdmf";
}

void writeXdmf(const std::string& xdmfPath, const std::string& h5Name,
               int numNodes, int numElements, const std::string& fieldName) {
  std::ofstream out(xdmfPath);
  if (!out) {
    throw std::runtime_error("femheat: cannot open XDMF file: " + xdmfPath);
  }
  out << "<?xml version=\"1.0\" ?>\n"
      << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
      << "<Xdmf Version=\"2.0\">\n"
      << "  <Domain>\n"
      << "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n"
      << "      <Topology TopologyType=\"Triangle\" NumberOfElements=\""
      << numElements << "\">\n"
      << "        <DataItem Dimensions=\"" << numElements
      << " 3\" NumberType=\"Int\" Format=\"HDF\">" << h5Name
      << ":/mesh/connectivity</DataItem>\n"
      << "      </Topology>\n"
      << "      <Geometry GeometryType=\"XY\">\n"
      << "        <DataItem Dimensions=\"" << numNodes
      << " 2\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << h5Name
      << ":/mesh/coordinates</DataItem>\n"
      << "      </Geometry>\n"
      << "      <Attribute Name=\"" << fieldName
      << "\" AttributeType=\"Scalar\" Center=\"Node\">\n"
      << "        <DataItem Dimensions=\"" << numNodes
      << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << h5Name
      << ":/solution/" << fieldName << "</DataItem>\n"
      << "      </Attribute>\n"
      << "    </Grid>\n"
      << "  </Domain>\n"
      << "</Xdmf>\n";
}

}  // namespace

void Hdf5Writer::writeTriMesh(const std::string& h5Path, const TriMesh& mesh,
                              const Eigen::VectorXd& field,
                              const std::string& fieldName) {
  const int numNodes = mesh.numNodes();
  const int numElements = mesh.numElements();
  if (field.size() != static_cast<Eigen::Index>(numNodes)) {
    throw std::runtime_error("femheat: field size does not match node count");
  }

  const auto& nodes = mesh.nodes();
  std::vector<std::vector<double>> coords(static_cast<std::size_t>(numNodes),
                                          std::vector<double>(2));
  for (int i = 0; i < numNodes; ++i) {
    const auto ui = static_cast<std::size_t>(i);
    coords[ui][0] = nodes[ui].x;
    coords[ui][1] = nodes[ui].y;
  }

  const auto& tris = mesh.triangles();
  std::vector<std::vector<int>> conn(static_cast<std::size_t>(numElements),
                                     std::vector<int>(3));
  for (int e = 0; e < numElements; ++e) {
    const auto& t = tris[static_cast<std::size_t>(e)];
    conn[static_cast<std::size_t>(e)] = {t[0], t[1], t[2]};
  }

  const std::vector<double> values(field.data(), field.data() + field.size());

  HighFive::File file(h5Path, HighFive::File::Overwrite);
  HighFive::Group meshGroup = file.createGroup("mesh");
  meshGroup.createDataSet("coordinates", coords);
  meshGroup.createDataSet("connectivity", conn);
  HighFive::Group solutionGroup = file.createGroup("solution");
  solutionGroup.createDataSet(fieldName, values);

  writeXdmf(xdmfPathFor(h5Path), h5BaseName(h5Path), numNodes, numElements,
            fieldName);
}

}  // namespace femheat
