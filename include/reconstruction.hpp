#pragma once

#include <cmath>
#include <vector>

#include "types.hpp"

namespace euler1d {

inline double limiterValue(Limiter limiter, double r) {
  switch (limiter) {
    case Limiter::None:
      return 1.0;
    case Limiter::Minmod:
      return std::max(0.0, std::min(1.0, r));
    case Limiter::VanLeer:
      return (r + std::fabs(r)) / (1.0 + std::fabs(r));
  }
  return 1.0;
}

inline void reconstructFaces(const std::vector<Conserved>& U,
                             std::vector<CellFaces>& Ufaces, Scheme scheme,
                             Limiter limiter) {
  const int n = static_cast<int>(U.size());
  constexpr double eps = 1e-12;

  if (scheme == Scheme::Constant) {
    for (int i = 0; i < n; ++i) {
      Ufaces[i][static_cast<int>(Face::West)] = U[i];
      Ufaces[i][static_cast<int>(Face::East)] = U[i];
    }
    return;
  }

  Ufaces[0][static_cast<int>(Face::West)] = U[0];
  Ufaces[0][static_cast<int>(Face::East)] = U[0];
  Ufaces[n - 1][static_cast<int>(Face::West)] = U[n - 1];
  Ufaces[n - 1][static_cast<int>(Face::East)] = U[n - 1];

  for (int i = 1; i < n - 1; ++i) {
    for (int v = 0; v < kNumVars; ++v) {
      const double dMinus = U[i][v] - U[i - 1][v];
      const double dPlus = U[i + 1][v] - U[i][v];
      const double r = dMinus / (dPlus + eps);
      const double phi = limiterValue(limiter, r);
      const double slope = phi * dPlus;

      Ufaces[i][static_cast<int>(Face::West)][v] = U[i][v] - 0.5 * slope;
      Ufaces[i][static_cast<int>(Face::East)][v] = U[i][v] + 0.5 * slope;
    }
  }
}

}  // namespace euler1d
