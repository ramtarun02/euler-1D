#pragma once

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "physics.hpp"
#include "types.hpp"

namespace euler1d {

inline std::string zeroPadded(int value) {
  std::ostringstream out;
  out << std::setfill('0') << std::setw(6) << value;
  return out.str();
}

inline void writeSolutionCsv(const std::vector<double>& x,
                             const std::vector<Conserved>& U,
                             const SimulationParams& simParams,
                             const std::string& outputPrefix = "solution") {
  const std::string fileName = outputPrefix + "_" + zeroPadded(simParams.numPoints) +
                               "_" + zeroPadded(simParams.timeStep) + ".csv";
  std::ofstream outputFile(fileName);
  outputFile << "x,rho,u,p\n";

  for (int i = 0; i < simParams.numPoints; ++i) {
    const double rho = U[i][0];
    const double u = U[i][1] / rho;
    const double p = pressureFromConserved(simParams.gamma, U[i]);
    outputFile << x[i] << ',' << rho << ',' << u << ',' << p << '\n';
  }
}

}  // namespace euler1d
