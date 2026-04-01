#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "types.hpp"

namespace euler1d {

inline double calcSpeedOfSound(double gamma, double rho, double p) {
  return std::sqrt(gamma * p / rho);
}

inline double pressureFromConserved(double gamma, const Conserved& U) {
  const double rho = U[0];
  const double u = U[1] / rho;
  const double kinetic = 0.5 * U[1] * u;
  return (gamma - 1.0) * (U[2] - kinetic);
}

inline Conserved fluxFromConserved(double gamma, const Conserved& U) {
  const double rho = U[0];
  const double u = U[1] / rho;
  const double p = pressureFromConserved(gamma, U);
  return {rho * u, p + rho * u * u, u * (U[2] + p)};
}

inline double maxWaveSpeed(double gamma, const Conserved& U) {
  const double rho = U[0];
  const double u = U[1] / rho;
  const double p = pressureFromConserved(gamma, U);
  return std::fabs(u) + calcSpeedOfSound(gamma, rho, p);
}

inline Conserved rusanovFlux(double gamma, const Conserved& left,
                            const Conserved& right) {
  const Conserved fluxL = fluxFromConserved(gamma, left);
  const Conserved fluxR = fluxFromConserved(gamma, right);
  const double sMax = std::max(maxWaveSpeed(gamma, left), maxWaveSpeed(gamma, right));

  Conserved F{};
  for (int j = 0; j < kNumVars; ++j) {
    F[j] = 0.5 * ((fluxL[j] + fluxR[j]) - sMax * (right[j] - left[j]));
  }
  return F;
}

}  // namespace euler1d
