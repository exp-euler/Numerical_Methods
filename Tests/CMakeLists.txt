# line below builds serial executable (remove mpi functions from main.cpp)
set(CPP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
)


# next lines build serial executables
if(CMAKE_BUILD_TYPE MATCHES Serial)
    set(CPP_SOURCE
        ${CPP_SOURCE}
        "${PROJECT_SOURCE_DIR}/Tests/matFunTest.cpp"
    )
    add_executable(matFunTest ${CPP_SOURCE})
endif(CMAKE_BUILD_TYPE MATCHES Serial)

# next lines build mpi-openmp executables
if(CMAKE_BUILD_TYPE MATCHES HybridCPU)
    set(CPP_SOURCE
        ${CPP_SOURCE}
        "${PROJECT_SOURCE_DIR}/Tests/matrixTest.cpp"
    )
    add_executable(matrixTest ${CPP_SOURCE})
    target_link_libraries(matrixTest PUBLIC ${MPI_LIBRARIES})
    target_link_libraries(matrixTest PUBLIC OpenMP::OpenMP_CXX)
endif(CMAKE_BUILD_TYPE MATCHES HybridCPU)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
