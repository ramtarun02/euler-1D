cmake_minimum_required(VERSION 3.10)
project(euler)

# ensure that the C++ 2017 standard is used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)

add_library(euler_core
	src/solver.cpp
)

target_include_directories(euler_core
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(euler_core
	PRIVATE -O3 -Wall -Wextra -Wpedantic
)

# Modular optimized solver (default executable)
add_executable(euler
	src/main.cpp
)
target_link_libraries(euler PRIVATE euler_core)

# Keep original monolithic file available unchanged
add_executable(euler_legacy euler.cpp)
