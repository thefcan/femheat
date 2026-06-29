#include "femheat/mesh1d.hpp"

#include "femheat/line2.hpp"

namespace femheat {

Mesh1D::Mesh1D(double a, double b, int numElements) {
  coords_.resize(static_cast<std::size_t>(numElements) + 1);
  const double h = (b - a) / numElements;
  for (int i = 0; i <= numElements; ++i) {
    coords_[static_cast<std::size_t>(i)] = a + i * h;
  }
}

std::vector<std::unique_ptr<IElement>> Mesh1D::buildElements() const {
  std::vector<std::unique_ptr<IElement>> elements;
  elements.reserve(static_cast<std::size_t>(numElements()));
  for (int e = 0; e < numElements(); ++e) {
    elements.push_back(std::make_unique<Line2>(
        e, e + 1, coords_[static_cast<std::size_t>(e)],
        coords_[static_cast<std::size_t>(e) + 1]));
  }
  return elements;
}

}  // namespace femheat
