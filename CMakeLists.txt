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

# Find dependencies
find_package(Eigen3 REQUIRED)

########################
## Library definition ##
########################
add_library(knn
  ${NABO_SRC}
)
target_include_directories(knn PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(knn PRIVATE
  Eigen3::Eigen
)
set_target_properties(knn PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(${PROJECT_NAME}::knn ALIAS knn)

#############
## Install ##
#############
install(
  TARGETS
    knn
  EXPORT ${PROJECT_NAME}_Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)
install(DIRECTORY include/
  DESTINATION include
)
install(FILES ${CMAKE_BINARY_DIR}/compile_commands.json
  DESTINATION .
  OPTIONAL
)
install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE ${PROJECT_NAME}Config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION share/${PROJECT_NAME}
)
install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE Find${PROJECT_NAME}.cmake
  DESTINATION share/${PROJECT_NAME}
)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION share/${PROJECT_NAME}
)