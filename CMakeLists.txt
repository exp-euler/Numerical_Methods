project(numerical_methods)
CMAKE_MINIMUM_REQUIRED(VERSION 3.21)

#########################################################################
#                               Compiler Setup
#########################################################################

# 1. 'O3': enable all optimizations
# 2. '-Wall': enable all warnings
SET(CMAKE_CXX_FLAGS "-std=c++17 -g -O3 -Wall")

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Serial CACHE STRING "Choose the type of build: Serial SharedCPU DistCPU HybridCPU GPU" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


#########################################################################
#                               Parallel Computing
#########################################################################

if(CMAKE_BUILD_TYPE MATCHES Serial)
    SET(CMAKE_CXX_COMPILER "g++")
    add_definitions(-DSERIAL)
endif(CMAKE_BUILD_TYPE MATCHES Serial)
if(CMAKE_BUILD_TYPE MATCHES SharedCPU)
    SET(CMAKE_CXX_COMPILER "g++")
    find_package(OpenMP REQUIRED)
    add_definitions(-DSHARED_CPU)
endif(CMAKE_BUILD_TYPE MATCHES SharedCPU)
if(CMAKE_BUILD_TYPE MATCHES DistCPU)
    SET(CMAKE_CXX_COMPILER "mpic++")
    find_package(MPI REQUIRED)
    add_definitions(-DDIST_CPU)
endif(CMAKE_BUILD_TYPE MATCHES DistCPU)
if(CMAKE_BUILD_TYPE MATCHES HybridCPU)
    SET(CMAKE_CXX_COMPILER "mpic++")
    find_package(MPI REQUIRED)
    find_package(OpenMP REQUIRED)
    add_definitions(-DHYBRID_CPU)
endif(CMAKE_BUILD_TYPE MATCHES HybridCPU)
if(CMAKE_BUILD_TYPE MATCHES GPU)
    SET(CMAKE_CXX_COMPILER "nvcc")
    add_definitions(-DGPU)
endif(CMAKE_BUILD_TYPE MATCHES GPU)


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
