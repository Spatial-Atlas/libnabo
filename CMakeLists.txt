cmake_minimum_required(VERSION 3.15)
project(nabo VERSION 1.0.0 LANGUAGES CXX)

# Source and header files.
set(NABO_SRC
  src/${PROJECT_NAME}/nabo.cpp
  src/${PROJECT_NAME}/kdtree_cpu.cpp
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20) # Use C++20 standard
set(CMAKE_CXX_STANDARD_REQUIRED ON) # Ensures that the specified C++ standard (set elsewhere) is strictly required for the project.
set(CMAKE_CXX_EXTENSIONS OFF) # Disables compiler-specific language extensions, enforcing strict compliance with the chosen C++ standard.

# Export compile commands for CI/CD tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define compiler options
add_compile_options(-Wall -Wextra -Wpedantic -Werror=return-type -Wunused -Wunused-function -Wunused-variable -Wunreachable-code)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-O3)
endif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")

# Build type options:
#	* catkin.
#	* cmake.
#	* ament (future).
set(BUILD_TYPE "cmake")
if (BUILD_TYPE STREQUAL "catkin")
  include(CatkinBuild.cmake)
else()
  # For a pure-cmake build, follow the instructions on the package.xml file, which
  # specify how to set the buildtool to cmake and export the package accordingly.
  include(PureCMakeBuild.cmake)
endif()
