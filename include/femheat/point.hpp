#ifndef FEMHEAT_POINT_HPP
#define FEMHEAT_POINT_HPP

namespace femheat {

/// A point in space. For 1D problems only `x` is used; `y` is there to support
/// the 2D triangle elements added in later milestones.
struct Point {
  double x{0.0};
  double y{0.0};
};

}  // namespace femheat

#endif  // FEMHEAT_POINT_HPP
