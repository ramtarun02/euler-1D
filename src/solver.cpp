#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "io.hpp"
#include "physics.hpp"
#include "reconstruction.hpp"

namespace euler1d {

namespace {

Conserved primitiveToConserved(double gamma, double rho, double u, double p) {
  return {rho, rho * u, p / (gamma - 1.0) + 0.5 * rho * u * u};
}

void applyBoundaryConditions(std::vector<Conserved>& U, const InitialCondition& ic) {
  const int n = static_cast<int>(U.size());

  U[0][0] = U[1][0];
  U[0][1] = ic.rhoLeft * ic.uLeft;
  U[0][2] = U[1][2];

  U[n - 1][0] = U[n - 2][0];
  U[n - 1][1] = ic.rhoRight * ic.uRight;
  U[n - 1][2] = U[n - 2][2];
}

}  // namespace

void runSimulation(SimulationParams simParams, const InitialCondition& ic,
                   Scheme scheme, Limiter limiter,
                   const std::string& outputPrefix) {
  simParams.dx = simParams.domainLength / (simParams.numPoints - 1);
  simParams.time = 0.0;
  simParams.timeStep = 0;

  const int n = simParams.numPoints;
  const double invDx = 1.0 / simParams.dx;

  std::vector<double> x(n);
  std::vector<Conserved> U(n);
  std::vector<CellFaces> Ufaces(n);
  std::vector<CellFaces> Ffaces(n);

  for (int i = 0; i < n; ++i) {
    x[i] = i * simParams.dx;
    if (x[i] < ic.xDiscontinuity) {
      U[i] = primitiveToConserved(simParams.gamma, ic.rhoLeft, ic.uLeft, ic.pLeft);
    } else {
      U[i] = primitiveToConserved(simParams.gamma, ic.rhoRight, ic.uRight, ic.pRight);
    }
  }

  while (simParams.time < simParams.endTime) {
    const std::vector<Conserved> UOld = U;

    double speedMax = 0.0;
    for (int i = 0; i < n; ++i) {
      speedMax = std::max(speedMax, maxWaveSpeed(simParams.gamma, U[i]));
    }

    double dt = (simParams.CFL * simParams.dx) / speedMax;
    dt = std::min(dt, simParams.endTime - simParams.time);

    reconstructFaces(U, Ufaces, scheme, limiter);

    for (int i = 1; i < n - 1; ++i) {
      Ffaces[i][static_cast<int>(Face::West)] =
          rusanovFlux(simParams.gamma,
                      Ufaces[i - 1][static_cast<int>(Face::East)],
                      Ufaces[i][static_cast<int>(Face::West)]);

      Ffaces[i][static_cast<int>(Face::East)] =
          rusanovFlux(simParams.gamma,
                      Ufaces[i][static_cast<int>(Face::East)],
                      Ufaces[i + 1][static_cast<int>(Face::West)]);
    }

    for (int i = 1; i < n - 1; ++i) {
      for (int j = 0; j < kNumVars; ++j) {
        const double dF = Ffaces[i][static_cast<int>(Face::East)][j] -
                          Ffaces[i][static_cast<int>(Face::West)][j];
        U[i][j] = UOld[i][j] - dt * invDx * dF;
      }
    }

    applyBoundaryConditions(U, ic);
    writeSolutionCsv(x, U, simParams, outputPrefix);

    std::cout << "Current time: " << std::scientific << std::setw(10)
              << std::setprecision(3) << simParams.time;
    std::cout << ", End time: " << std::scientific << std::setw(10)
              << std::setprecision(3) << simParams.endTime;
    std::cout << ", Current time step: " << std::fixed << std::setw(7)
              << simParams.timeStep << "\r";

    simParams.time += dt;
    simParams.timeStep++;
  }

  std::cout << "\nSimulation finished" << std::endl;
}

}  // namespace euler1d
