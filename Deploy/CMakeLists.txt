# next lines build serial executables
if(CMAKE_BUILD_TYPE MATCHES Serial)
    add_executable(scalarIVP ${CPP_SOURCE} ./scalarIVP.cpp)
endif(CMAKE_BUILD_TYPE MATCHES Serial)

# next lines build mpi-openmp executables
if(CMAKE_BUILD_TYPE MATCHES HybridCPU)
    add_executable(scalarIVP ${CPP_SOURCE} ./scalarIVP.cpp)
    target_link_libraries(scalarIVP PUBLIC ${MPI_LIBRARIES})
    target_link_libraries(scalarIVP PUBLIC OpenMP::OpenMP_CXX)
endif(CMAKE_BUILD_TYPE MATCHES HybridCPU)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
