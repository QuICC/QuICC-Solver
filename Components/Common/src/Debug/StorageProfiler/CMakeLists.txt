if(QUICC_STORAGEPROFILE)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    StorageProfiler.cpp
    StorageProfilerTools.cpp
    )
endif(QUICC_STORAGEPROFILE)
