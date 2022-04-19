if(QUICC_PROFILE)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    BreakPoint.cpp
    ProfilerBase.cpp
    Tools.cpp
    )

  # Add in the parallel or serial profiler depending on setup
  if(QUICC_MPI)
    quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
      MpiProfiler.cpp
      )
  else(QUICC_MPI)
    quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
      SerialProfiler.cpp
      )
  endif(QUICC_MPI)
endif(QUICC_PROFILE)
