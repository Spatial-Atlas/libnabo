# Find catkin macros and libraries
find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)

# Catkin package macro
catkin_package(
  INCLUDE_DIRS
    ${PROJECT_NAME}
    ${CMAKE_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
  LIBRARIES
    ${PROJECT_NAME}
)

########################
## Library definition ##
########################
# Nabo
add_library(${PROJECT_NAME}
  ${NABO_SRC}
)

target_include_directories(${PROJECT_NAME}
  PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/nabo
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  Eigen3::Eigen
)

#############
## Install ##
#############
install(
  TARGETS
    nabo
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(
  DIRECTORY
    ${CMAKE_SOURCE_DIR}/nabo/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}/nabo
)