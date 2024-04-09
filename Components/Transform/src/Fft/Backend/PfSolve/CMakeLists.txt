if(QUICC_FFT_WORLAND STREQUAL "PfSolve")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IWorlandBackend.cpp
    WorlandProjector.cpp
    WorlandIntegrator.cpp
  )
  set(_use_PFSOLVE ON)
endif()

if(_use_PFSOLVE)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IPfSolveBackend.cpp
    Tools.cpp
  )
endif()
