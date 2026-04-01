#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

enum SCHEME { CONSTANT = 0, MUSCL };
enum LIMITER { NONE = 0, MINMOD, VANLEER };
enum FACE { WEST = 0, EAST };

struct simulationParams {
  int numPoints;
  double gamma;
  double domainLength;
  double CFL;
  double endTime;
  double dx;
  double time;
  int timeStep;
};

double calcSpeedOfSound(double gamma, double rho, double p) {
  return std::sqrt(gamma * p / rho);
}

int main() {
  // ------------------Prep-Processing Block -------------------
  // 1. Setting up the case simParams
  auto numericalScheme = SCHEME::MUSCL;
  auto limiter = LIMITER::VANLEER;

  simulationParams simParams;

  simParams.domainLength = 1.0;
  simParams.CFL = 0.10;
  simParams.gamma = 1.4;
  simParams.numPoints = 101;

  simParams.dx = simParams.domainLength / (simParams.numPoints - 1);
  simParams.endTime = 0.2;
  simParams.time = 0.0;
  simParams.timeStep = 0;

  // 2. Allocating Memory for Data Structures
  std::vector<double> x(simParams.numPoints);

  std::vector<std::array<double, 3>> U(simParams.numPoints);

  // We use a std::vector to have an entry available for each cell in the
  // domain. Then, we want to store, for each cell, an std::array containing 3
  // elements.
  std::vector<std::array<std::array<double, 3>, 2>> Ufaces(simParams.numPoints);
  std::vector<std::array<std::array<double, 3>, 2>> Ffaces(simParams.numPoints);

  // Create/Read Mesh
  // Fill the x vector once the memory is allocated

  for (int i = 0; i < simParams.numPoints; i++) {
    x[i] = i * simParams.dx;
  }

  // Initialise the Solution to our cells
  //
  //
  double u = 0.0;
  double rho = 0.0;
  double p = 0.0;

  for (int i = 0; i < simParams.numPoints; i++) {
    if (x[i] < 0.5) {
      u = 0.0;
      rho = 1.0;
      p = 1.0;

    } else {
      rho = 0.125;
      u = 0;
      p = 0.1;
    }

    U[i][0] = rho;
    U[i][1] = rho * u;
    U[i][2] = p / (simParams.gamma - 1) + 0.5 * rho * std::pow(u, 2);
  }

  // -------------------------------Solver-------------------------
  while (simParams.time < simParams.endTime) {
    // Preparing solution update (store old solution and calculate stable
    auto UOld = U;

    // calculate stable time step
    double speedMax = 0.0;
    for (int i = 0; i < simParams.numPoints; i++) {
      // we need the primitive variables first to compute the wave speed (
      // based on speed of sound and local velocity)
      auto rho = U[i][0];
      auto u = U[i][1] / rho;
      auto p = (simParams.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
      // calculate wave speed for each cell
      double speedOfSound = calcSpeedOfSound(simParams.gamma, rho, p);
      if (speedOfSound + std::fabs(u) > speedMax)
        speedMax = speedOfSound + std::fabs(u);
    }
    double dt = (simParams.CFL * simParams.dx) / speedMax;

    if (numericalScheme == SCHEME::CONSTANT) {
      for (int i = 0; i < simParams.numPoints; i++) {
        for (int j = 0; j < 3; ++j) {
          Ufaces[i][FACE::WEST][j] = U[i][j];
          Ufaces[i][FACE::EAST][j] = U[i][j];
        }
      }
    } else if (numericalScheme == SCHEME::MUSCL) {
      // Use Lower order schemes near the boundaries
      for (int variable = 0; variable < 3; ++variable) {
        Ufaces[0][FACE::WEST][variable] = U[0][variable];
        Ufaces[0][FACE::EAST][variable] = U[0][variable];

        Ufaces[simParams.numPoints - 1][FACE::WEST][variable] =
            U[simParams.numPoints - 1][variable];
        Ufaces[simParams.numPoints - 1][FACE::EAST][variable] =
            U[simParams.numPoints - 1][variable];
      }

      // Using high resolutions MUSCL scheme on interior nodes
      for (int i = 1; i < simParams.numPoints - 1; i++) {
        for (int variable = 0; variable < 3; ++variable) {
          auto du_i_plus_half = U[i + 1][variable] - U[i][variable];
          auto du_i_minus_half = U[i][variable] - U[i - 1][variable];

          double rL = du_i_minus_half / (du_i_plus_half + 1e-8);
          double rR = du_i_plus_half / (du_i_minus_half + 1e-8);

          double psiL = 1.0;
          double psiR = 1.0;

          // apply limiters to make scheme TVD
          if (limiter == LIMITER::MINMOD) {
            psiL = std::max(0.0, std::min(1.0, rL));
            psiR = std::max(0.0, std::min(1.0, rR));
          } else if (limiter == LIMITER::VANLEER) {
            psiL = (rL + std::fabs(rL)) / (1.0 + std::fabs(rL));
            psiR = (rR + std::fabs(rR)) / (1.0 + std::fabs(rR));
          }

          Ufaces[i][FACE::WEST][variable] =
              U[i][variable] - 0.5 * psiL * du_i_plus_half;
          Ufaces[i][FACE::EAST][variable] =
              U[i][variable] + 0.5 * psiR * du_i_minus_half;
        }
      }
    }

    // Compute Flux at the Faces

    std::array<double, 3> fluxL, fluxR;

    for (int i = 1; i < simParams.numPoints - 1; i++) {
      for (int face = FACE::WEST; face <= FACE::EAST; ++face) {
        int indexOffset = 0;
        if (face == FACE::WEST)
          indexOffset = 0;
        else if (face == FACE::EAST)
          indexOffset = 1;

        auto rhoL = Ufaces[i - 1 + indexOffset][FACE::EAST][0];
        auto uL = Ufaces[i - 1 + indexOffset][FACE::EAST][1] / rhoL;
        auto EL = Ufaces[i - 1 + indexOffset][FACE::EAST][2];
        auto pL = (simParams.gamma - 1) * (EL - 0.5 * rhoL * std::pow(uL, 2));
        auto aL = calcSpeedOfSound(simParams.gamma, rhoL, pL);

        auto rhoR = Ufaces[i + indexOffset][FACE::WEST][0];
        auto uR = Ufaces[i + indexOffset][FACE::WEST][1] / rhoR;
        auto ER = Ufaces[i + indexOffset][FACE::WEST][2];
        auto pR = (simParams.gamma - 1) * (ER - 0.5 * rhoR * std::pow(uR, 2));
        auto aR = calcSpeedOfSound(simParams.gamma, rhoR, pR);

        fluxL[0] = rhoL * uL;
        fluxL[1] = pL + rhoL * std::pow(uL, 2);
        fluxL[2] = uL * (EL + pL);

        fluxR[0] = rhoR * uR;
        fluxR[1] = pR + rhoR * std::pow(uR, 2);
        fluxR[2] = uR * (ER + pR);

        // Rusanov Riemann Solver (Local Lax Friedrichs)
        auto speedMax = std::max(std::fabs(uL) + aL, std::fabs(uR) + aR);
        for (int j = 0; j < 3; ++j) {
          const auto &qL = Ufaces[i - 1 + indexOffset][FACE::EAST][j];
          const auto &qR = Ufaces[i + indexOffset][FACE::WEST][j];
          const auto &fL = fluxL[j];
          const auto &fR = fluxR[j];

          Ffaces[i][face][j] = 0.5 * ((fL + fR) - speedMax * (qR - qL));
        }
      }
    }
    // Update Solution in Time
    for (int i = 1; i < simParams.numPoints - 1; i++) {
      for (int j = 0; j < 3; j++) {
        const auto &dF = Ffaces[i][FACE::EAST][j] - Ffaces[i][FACE::WEST][j];
        U[i][j] = UOld[i][j] - (dt / simParams.dx) * dF;
      }
    }

    auto rhoL = 1.0;
    auto uL = 0.0;

    auto rhoR = 0.125;
    auto uR = 0.0;

    U[0][0] = U[1][0];
    U[0][1] = rhoL * uL;
    U[0][2] = U[1][2];
    U[simParams.numPoints - 1][0] = U[simParams.numPoints - 2][0];
    U[simParams.numPoints - 1][1] = rhoR * uR;
    U[simParams.numPoints - 1][2] = U[simParams.numPoints - 2][2];

    // output solution to a csv file for plotting
    std::ofstream outputFile;
    // Convert timestep and points into a 6 digit string with leading zeros
    std::ostringstream timeStepTemp, pointsTemp;

    timeStepTemp << std::setfill('0') << std::setw(6);
    timeStepTemp << simParams.timeStep;
    auto timeStep = timeStepTemp.str();

    pointsTemp << std::setfill('0') << std::setw(6);
    pointsTemp << simParams.numPoints;
    auto points = pointsTemp.str();

    outputFile.open("solution_" + points + "_" + timeStep + ".csv");
    outputFile << "x,rho,u,p" << std::endl;
    for (int i = 0; i < simParams.numPoints; i++) {
      auto rho = U[i][0];
      auto u = U[i][1] / rho;
      auto p = (simParams.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
      outputFile << x[i] << "," << rho << "," << u << "," << p << std::endl;
    }
    outputFile.close();

    // output current time step information to screen
    std::cout << "Current time: " << std::scientific << std::setw(10)
              << std::setprecision(3) << simParams.time;

    std::cout << ", End time: " << std::scientific << std::setw(10)
              << std::setprecision(3) << simParams.endTime;
    std::cout << ", Current time step: " << std::fixed << std::setw(7)
              << simParams.timeStep;
    std::cout << "\r";

    // Increment solution time and time step
    simParams.time += dt;
    simParams.timeStep++;
  }

  std::cout << "\nSimulation finished" << std::endl;

  return 0;
}
