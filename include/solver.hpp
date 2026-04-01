#pragma once

#include <string>

#include "types.hpp"

namespace euler1d {

void runSimulation(SimulationParams simParams, const InitialCondition& ic,
                   Scheme scheme, Limiter limiter,
                   const std::string& outputPrefix = "solution");

}  // namespace euler1d
