# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT(polar2Cartesian)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})


ADD_EXECUTABLE(polar2Cartesian Polar2Cartesian.cxx )

TARGET_LINK_LIBRARIES(polar2Cartesian  ${ITK_LIBRARIES})
