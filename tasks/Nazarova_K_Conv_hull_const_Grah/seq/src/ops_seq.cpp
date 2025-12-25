#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <vector>

#include "Nazarova_K_Conv_hull_const_Grah/common/include/common.hpp"
#include "Nazarova_K_Conv_hull_const_Grah/seq/include/ops_seq.hpp"

namespace nazarova_k_conv_hull_const_grah_processes {
namespace {

inline bool LessPivot(const Point& a, const Point& b) {
  if (a.y != b.y) {
    return a.y < b.y;
  }
  return a.x < b.x;
}

inline bool LessXY(const Point& a, const Point& b) {
  if (a.x != b.x) {
    return a.x < b.x;
  }
  return a.y < b.y;
}

inline std::int64_t Cross(const Point& o, const Point& a, const Point& b) {
  const std::int64_t ax = static_cast<std::int64_t>(a.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t ay = static_cast<std::int64_t>(a.y) - static_cast<std::int64_t>(o.y);
  const std::int64_t bx = static_cast<std::int64_t>(b.x) - static_cast<std::int64_t>(o.x);
  const std::int64_t by = static_cast<std::int64_t>(b.y) - static_cast<std::int64_t>(o.y);
  return (ax * by) - (ay * bx);
}

inline std::int64_t Dist2(const Point& a, const Point& b) {
  const std::int64_t dx = static_cast<std::int64_t>(a.x) - static_cast<std::int64_t>(b.x);
  const std::int64_t dy = static_cast<std::int64_t>(a.y) - static_cast<std::int64_t>(b.y);
  return (dx * dx) + (dy * dy);
}

std::vector<Point> GrahamScan(std::vector<Point> pts) {
  std::ranges::sort(pts, LessXY);
  pts.erase(std::ranges::unique(pts).begin(), pts.end());

  if (pts.size() <= 1U) {
    return pts;
  }

  // Pivot: lowest (y, x)
  const auto pivot_it = std::ranges::min_element(pts, LessPivot);
  std::iter_swap(pts.begin(), pivot_it);
  const Point pivot = pts.front();

  auto angle_less = [&pivot](const Point& a, const Point& b) {
    const std::int64_t cr = Cross(pivot, a, b);
    if (cr != 0) {
      return cr > 0;  // CCW first
    }
    // Same ray: keep closer first, so farther one survives when popping collinear points.
    return Dist2(pivot, a) < Dist2(pivot, b);
  };

  std::ranges::sort(std::next(pts.begin()), pts.end(), angle_less);

  if (pts.size() == 2U) {
    return pts;
  }

  std::vector<Point> hull;
  hull.reserve(pts.size());
  hull.push_back(pts[0]);
  hull.push_back(pts[1]);

  for (std::size_t i = 2; i < pts.size(); i++) {
    while (hull.size() >= 2U &&
           Cross(hull[hull.size() - 2U], hull[hull.size() - 1U], pts[i]) <= static_cast<std::int64_t>(0)) {
      hull.pop_back();
    }
    hull.push_back(pts[i]);
  }

  // Degenerate case: all points collinear -> keep endpoints
  if (hull.size() == 1U) {
    return hull;
  }
  if (hull.size() == 2U) {
    if (LessPivot(hull[1], hull[0])) {
      std::swap(hull[0], hull[1]);
    }
  }

  return hull;
}

}  // namespace

NazarovaKConvHullConstGrahSEQ::NazarovaKConvHullConstGrahSEQ(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool NazarovaKConvHullConstGrahSEQ::ValidationImpl() {
  // Any number of points is allowed (including empty).
  return true;
}

bool NazarovaKConvHullConstGrahSEQ::PreProcessingImpl() {
  points_ = GetInput().points;
  GetOutput().clear();
  return true;
}

bool NazarovaKConvHullConstGrahSEQ::RunImpl() {
  GetOutput() = GrahamScan(points_);
  return true;
}

bool NazarovaKConvHullConstGrahSEQ::PostProcessingImpl() {
  // Hull cannot have more points than input (after duplicate removal).
  return GetOutput().size() <= GetInput().points.size();
}

}  // namespace nazarova_k_conv_hull_const_grah_processes
