# line below builds serial executable (remove mpi functions from main.cpp)
set(CPP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/TimeIntegration.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
    "${PROJECT_SOURCE_DIR}/src/Matrix.cpp"
    "${PROJECT_SOURCE_DIR}/src/Vector.cpp"
    "${PROJECT_SOURCE_DIR}/Tests/matrixTest.cpp"
)

# next two lines build mpi executable
# add_executable(main ${CPP_SOURCE})
add_executable(matrixTest ${CPP_SOURCE})
target_link_libraries(matrixTest PUBLIC ${MPI_LIBRARIES})
target_link_libraries(matrixTest PUBLIC OpenMP::OpenMP_CXX)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
