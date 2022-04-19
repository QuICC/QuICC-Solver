quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Precision.cpp
  Profiler.cpp
  )

if(QUICC_LARGEIO STREQUAL "HDF5")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    Hdf5.cpp
    )
endif(QUICC_LARGEIO STREQUAL "HDF5")

if(QUICC_MPI)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    Mpi.cpp
    )
endif(QUICC_MPI)
