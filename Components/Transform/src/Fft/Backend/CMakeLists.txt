if(QUICC_FFT_MIXED STREQUAL "cuFFT")
  set(_useCUDA ON)
else()
  set(_useCUDA OFF)
endif()

quicc_target_cuda_sources(_useCUDA ${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  MixedProjector.cpp
  MixedIntegrator.cpp
)

if(QUICC_FFT_COMPLEX STREQUAL "cuFFT")
  set(_useCUDA ON)
else()
  set(_useCUDA OFF)
endif()

quicc_target_cuda_sources(_useCUDA ${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ComplexProjector.cpp
  ComplexIntegrator.cpp
)

if(QUICC_FFT_CHEBYSHEV STREQUAL "cuFFT")
  set(_useCUDA ON)
else()
  set(_useCUDA OFF)
endif()

quicc_target_cuda_sources(_useCUDA ${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ChebyshevProjector.cpp
  ChebyshevIntegrator.cpp
  ChebyshevEnergy.cpp
)

if(QUICC_FFT_WORLAND STREQUAL "cuFFT")
  set(_useCUDA ON)
else()
  set(_useCUDA OFF)
endif()

quicc_target_cuda_sources(_useCUDA ${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  WorlandProjector.cpp
  WorlandIntegrator.cpp
)

add_subdirectory(Fftw)
add_subdirectory(PfSolve)
add_subdirectory(CuFft)
