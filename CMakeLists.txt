
cmake_minimum_required(VERSION 3.20)

project(vascpp)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)
find_package(ParMETIS REQUIRED)

SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${PARMETIS_INCLUDE_DIR})

add_subdirectory(libs)
add_subdirectory(app)
