#include "femheat/assembler.hpp"

namespace femheat {

LinearSystem Assembler::assemble(int numNodes,
                                 const std::vector<std::unique_ptr<IElement>>& elements,
                                 const Material& material,
                                 const std::function<double(const Point&)>& source) {
  LinearSystem system;
  system.f = Eigen::VectorXd::Zero(numNodes);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(elements.size() * 9);  // up to 3x3 per element (2D)

  for (const auto& element : elements) {
    const Eigen::MatrixXd ke = element->elementStiffness(material);
    const Eigen::VectorXd fe = element->elementLoad(source(element->centroid()));
    const std::vector<int> dofs = element->dofs();

    for (std::size_t a = 0; a < dofs.size(); ++a) {
      const auto ia = static_cast<Eigen::Index>(a);
      system.f(dofs[a]) += fe(ia);
      for (std::size_t b = 0; b < dofs.size(); ++b) {
        triplets.emplace_back(dofs[a], dofs[b], ke(ia, static_cast<Eigen::Index>(b)));
      }
    }
  }

  system.K.resize(numNodes, numNodes);
  // setFromTriplets sums duplicate (i, j) contributions — exactly the scatter
  // step of finite element assembly.
  system.K.setFromTriplets(triplets.begin(), triplets.end());
  return system;
}

}  // namespace femheat
