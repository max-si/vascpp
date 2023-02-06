# set CMake version
cmake_minimum_required(VERSION 3.20)

# project name
project(vascpp)

# Find packages
find_package(HDF5 REQUIRED)
#add_subdirectory(HighFive)
find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)
# Echo trilinos build info just for fun
# MESSAGE("\nFound Trilinos!  Here are the details: ")
# MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
# MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
# MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
# MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
# MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
# MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
# MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
# MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
# MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
# MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
# MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
# MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
# MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
# MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
# MESSAGE("End of Trilinos details\n")

# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# specify directories
include_directories(${PROJECT_SOURCE_DIR})
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
#LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

# set sources
set(SOURCE_EXE 
	main.cpp coordinate.cpp network.cpp vessel.cpp 
	RootVesselExporter.cpp SubNetworkExporter.cpp GeneratorSharedHdfFunctions.cpp 
	ImportExportCommon.cpp VesselNetworkExporter_HDF5.cpp 
	PressureListCreator.cpp Triplet.cpp EpetraTutorials.cpp)

# executable
add_executable(vascpp ${SOURCE_EXE})

# link libraries
target_link_libraries (vascpp PUBLIC HDF5::HDF5 MPI::MPI_CXX)
target_link_libraries(vascpp PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
