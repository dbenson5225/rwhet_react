cmake_minimum_required(VERSION 2.8.12)
project(rwhet_react Fortran)

# CMake Configuration variables
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the build type: Debug Release" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Setup dependencies
set(PhreeqcRM_ROOT_DIR /usr/local/share/doc/phreeqcrm/src)
find_package(PhreeqcRM REQUIRED)
include_directories(${PhreeqcRM_INCLUDE_DIR})

message("PhreeqcRM_LIBRARIES ${PhreeqcRM_LIBRARIES}")

add_executable(tlp_exe test_lapack.f90)

target_link_libraries(tlp_exe ${LAPACK_LIBRARIES})
