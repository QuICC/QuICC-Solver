if(QUICC_USE_MPI)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    MpiFramework.cpp
    )
endif()
