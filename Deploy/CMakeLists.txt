# line below builds serial executable (remove mpi functions from main.cpp)
set(CPP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/Tableaus.cpp"
    #"${PROJECT_SOURCE_DIR}/src/Matrix.cpp"
    #"${PROJECT_SOURCE_DIR}/src/Vector.cpp"
    "${PROJECT_SOURCE_DIR}/Deploy/scalarIVP.cpp"
)
add_executable(scalarIVP ${CPP_SOURCE})

# next two lines build mpi executable
# add_executable(main ${CPP_SOURCE})
#target_link_libraries(main PUBLIC ${MPI_LIBRARIES})
#target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)

# checking where is mpi.h
#message(STATUS "mpi include dir: ${MPI_INCLUDE_PATH}")
