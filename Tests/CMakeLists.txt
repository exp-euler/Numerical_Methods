# next lines build serial executables
if(CMAKE_BUILD_TYPE MATCHES Serial)
    add_executable(matFunTest ${CPP_SOURCE} ./matFunTest.cpp)
    add_executable(vectorIVP ${CPP_SOURCE} ./vectorIVP.cpp)
endif(CMAKE_BUILD_TYPE MATCHES Serial)

# next lines build mpi-openmp executables
if(CMAKE_BUILD_TYPE MATCHES HybridCPU)
    add_executable(matrixTest ${CPP_SOURCE} ./matrixTest.cpp)
    target_link_libraries(matrixTest PUBLIC ${MPI_LIBRARIES})
    target_link_libraries(matrixTest PUBLIC OpenMP::OpenMP_CXX)
endif(CMAKE_BUILD_TYPE MATCHES HybridCPU)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
