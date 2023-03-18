project(numerical_methods)
cmake_minimum_required(VERSION 3.21)

#########################################################################
#                               Compiler Setup
#########################################################################

# use the C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# any compiler flags that should be used
#
# 1. 'O3': enable all optimizations
# 2. '-Wall': enable all warnings
set(CMAKE_CXX_FLAGS "-O3 -Wall")

#########################################################################
#                               Parallel Computing
#########################################################################

# parallel programing
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

################################################################################

# specify the non-header source files that should be compiled
# use GLOB-ing and * to select many files at once
# PROJECT_SOURCE_DIR is a built-in having the current dir of the project

#file(GLOB CPP_SOURCE
#    "${PROJECT_SOURCE_DIR}/src/*.cpp" #i.e. 'src/main.cpp', but not 'src/d/f.cpp'
#    "${PROJECT_SOURCE_DIR}/src/*/*.cpp" #i.e. 'src/d/f.cpp'
#)

################################################################################

#########################################################################
#                          Executables Directory
#########################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#########################################################################
#                       Source Files and Header Files
#########################################################################

#set(CPP_SOURCE
#    "${PROJECT_SOURCE_DIR}/src/TimeIntegration.cpp"
#    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
#    "${PROJECT_SOURCE_DIR}/src/Matrix.cpp"
#    "${PROJECT_SOURCE_DIR}/src/Vector.cpp"
#    "${PROJECT_SOURCE_DIR}/src/main.cpp"
#)

# where to find the header files
include_directories(${PROJECT_SOURCE_DIR}/src
                    ${MPI_INCLUDE_PATH})


#add_executable(main ${CPP_SOURCE})
#target_link_libraries(main PUBLIC ${MPI_LIBRARIES})
#target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
#########################################################################
#                Main files that use the functionalities
#########################################################################
include_directories(${PROJECT_SOURCE_DIR}/Deploy)
add_subdirectory(${PROJECT_SOURCE_DIR}/Deploy)

#########################################################################
#                Main files that test the functionalities
#########################################################################
include_directories(${PROJECT_SOURCE_DIR}/Tests)
add_subdirectory(${PROJECT_SOURCE_DIR}/Tests)
