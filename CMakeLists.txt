# Useful commands
# make clean will clean all intermediate files to re-compile
# make -j will compile in parallel


# suppress warnings
cmake_minimum_required(VERSION 3.21)

#
project(cpp_with_cmake)

# use the C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# any compiler flags that should be used
#
# 1. 'O3': enable all optimizations
# 2. '-Wall': enable all warnings
set(CMAKE_CXX_FLAGS "-O3 -Wall")

# specify the non-header source files that should be compiled
# use GLOB-ing and * to select many files at once
# PROJECT_SOURCE_DIR is a built-in having the current dir of the project

#file(GLOB CPP_SOURCE
#    "${PROJECT_SOURCE_DIR}/src/*.cpp" #i.e. 'src/main.cpp', but not 'src/d/f.cpp'
#    "${PROJECT_SOURCE_DIR}/src/*/*.cpp" #i.e. 'src/d/f.cpp'
#)

set(CPP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/TimeIntegration.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
    "${PROJECT_SOURCE_DIR}/src/main.cpp")

# where to find the header files
include_directories(${PROJECT_SOURCE_DIR}/src)

# name of the executable to be created, and which source files to use
add_executable(main ${CPP_SOURCE})
