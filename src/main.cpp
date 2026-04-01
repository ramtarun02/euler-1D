#include "solver.hpp"

int main() {
  euler1d::SimulationParams simParams;
  simParams.domainLength = 1.0;
  simParams.CFL = 0.10;
  simParams.gamma = 1.4;
  simParams.numPoints = 101;
  simParams.endTime = 0.2;

  euler1d::InitialCondition ic;

  euler1d::runSimulation(simParams, ic, euler1d::Scheme::Muscl,
                         euler1d::Limiter::VanLeer);
  return 0;
}
