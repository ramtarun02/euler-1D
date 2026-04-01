#pragma once

#include <array>

namespace euler1d {

enum class Scheme { Constant = 0, Muscl };
enum class Limiter { None = 0, Minmod, VanLeer };
enum class Face { West = 0, East = 1 };

constexpr int kNumVars = 3;

using Conserved = std::array<double, kNumVars>;
using CellFaces = std::array<Conserved, 2>;

struct SimulationParams {
  int numPoints = 101;
  double gamma = 1.4;
  double domainLength = 1.0;
  double CFL = 0.10;
  double endTime = 0.2;
  double dx = 0.0;
  double time = 0.0;
  int timeStep = 0;
};

struct InitialCondition {
  double xDiscontinuity = 0.5;
  double rhoLeft = 1.0;
  double uLeft = 0.0;
  double pLeft = 1.0;
  double rhoRight = 0.125;
  double uRight = 0.0;
  double pRight = 0.1;
};

} // namespace euler1d
