
add_executable(vessel_gen main.cpp)

#target_link_libraries(vessel_gen PRIVATE VascppLibrary ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES})

target_link_libraries(vessel_gen PRIVATE VascppLibrary ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} 
		${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${PARMETIS_LIBRARY} ${METIS_LIBRARY}
		Eigen3::Eigen) 

#${GSL_LIBRARIES} 