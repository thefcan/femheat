#ifndef FEMHEAT_MATERIAL_HPP
#define FEMHEAT_MATERIAL_HPP

namespace femheat {

/// Isotropic material holding the thermal conductivity k. A single global value
/// is stored for now; per-element conductivity can be layered on later.
class Material {
 public:
  explicit Material(double conductivity) : conductivity_(conductivity) {}

  double conductivity() const { return conductivity_; }

 private:
  double conductivity_;
};

}  // namespace femheat

#endif  // FEMHEAT_MATERIAL_HPP
