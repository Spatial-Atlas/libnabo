# Find catkin macros and libraries
find_package(Eigen3 REQUIRED)

########################
## Library definition ##
########################
add_library(knn
  ${NABO_SRC}
)
target_include_directories(knn
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(knn
  PUBLIC
    Eigen3::Eigen
)
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
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
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