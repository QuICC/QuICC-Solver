if(QUICC_FFT_MIXED STREQUAL "FFTW")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IMixedBackend.cpp
    MixedProjector.cpp
    MixedIntegrator.cpp
    )
  set(_use_FFTW ON)
endif()

if(QUICC_FFT_COMPLEX STREQUAL "FFTW")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IComplexBackend.cpp
    ComplexProjector.cpp
    ComplexIntegrator.cpp
  )
set(_use_FFTW ON)
endif()

if(QUICC_FFT_CHEBYSHEV STREQUAL "FFTW")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IChebyshevBackend.cpp
    ChebyshevProjector.cpp
    ChebyshevIntegrator.cpp
    ChebyshevEnergy.cpp
    DifferentialSolver.cpp
  )
set(_use_FFTW ON)
endif()

if(QUICC_FFT_WORLAND STREQUAL "FFTW")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IWorlandBackend.cpp
    WorlandProjector.cpp
    WorlandIntegrator.cpp
  )
set(_use_FFTW ON)
endif()

if(_use_FFTW)
quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IFftwBackend.cpp
  Library.cpp
  Tools.cpp
)
endif()
