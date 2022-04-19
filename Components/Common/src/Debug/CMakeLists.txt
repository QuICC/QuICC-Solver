if(QUICC_DEBUG)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    SerialDebugger.cpp
    )
endif(QUICC_DEBUG)

add_subdirectory(Profiler)
add_subdirectory(StorageProfiler)
