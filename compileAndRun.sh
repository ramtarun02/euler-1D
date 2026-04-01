#!/bin/bash

# Prepare the compilation (ensure that we have a clean build/ directory)
rm -rf build
mkdir build

# Compile the solver into an object file
g++ -std=c++17 -I. -c -O2 euler.cpp -o build/euler.o

# Link object file into an executable
g++ build/euler.o -o build/euler

# Run the solver
cd build
./euler
cd ..