# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(BHASagCalculator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Source files
set(SOURCES
    src/calculators/fea_core/fea_solver.cpp
    src/calculators/fea_core/trajectory.cpp
    src/calculators/fea_core/fea_results.cpp
    src/calculators/beam_elements/beam_element.cpp
    src/utils/file_readers.cpp
    src/calculators/main.cpp
)

# Executable
add_executable(bha_sag_calculator ${SOURCES})

# Copy data files to build directory
file(COPY src/calculators/data DESTINATION ${CMAKE_BINARY_DIR})