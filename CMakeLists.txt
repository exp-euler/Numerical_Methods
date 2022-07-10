# Useful commands
# make clean will clean all intermediate files to re-compile
# make -j will compile in parallel


# suppress warnings
cmake_minimum_required(VERSION 3.21)

project(numerical_methods)

# use the C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# parallel programing
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# any compiler flags that should be used
#
# 1. 'O3': enable all optimizations
# 2. '-Wall': enable all warnings
set(CMAKE_CXX_FLAGS "-O3 -Wall")

################################################################################

# specify the non-header source files that should be compiled
# use GLOB-ing and * to select many files at once
# PROJECT_SOURCE_DIR is a built-in having the current dir of the project

#file(GLOB CPP_SOURCE
#    "${PROJECT_SOURCE_DIR}/src/*.cpp" #i.e. 'src/main.cpp', but not 'src/d/f.cpp'
#    "${PROJECT_SOURCE_DIR}/src/*/*.cpp" #i.e. 'src/d/f.cpp'
#)

################################################################################

set(CPP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/TimeIntegration.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
    "${PROJECT_SOURCE_DIR}/src/Matrix.cpp"
    "${PROJECT_SOURCE_DIR}/src/Vector.cpp"
    "${PROJECT_SOURCE_DIR}/src/main.cpp"
)

# where to find the header files
include_directories(${PROJECT_SOURCE_DIR}/src
                    ${MPI_INCLUDE_PATH})

# line below builds serial executable (remove mpi functions from main.cpp)
#add_executable(main ${CPP_SOURCE})

# next two lines build mpi executable
add_executable(main ${CPP_SOURCE})
target_link_libraries(main PUBLIC ${MPI_LIBRARIES})
target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
