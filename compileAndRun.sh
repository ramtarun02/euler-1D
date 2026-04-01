#!/bin/bash

# Prepare the compilation (ensure that we have a clean build/ directory)
rm -rf build
mkdir build

# Configure and compile modular solver
cmake -S . -B build
cmake --build build --target euler -j

# Run the solver
(cd build && ./euler)
