if(QUICC_FFT_MIXED STREQUAL "cuFFT")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IMixedBackend.cu
    MixedProjector.cu
    MixedIntegrator.cu
    )
  set(_use_cuFFT ON)
endif()

if(QUICC_FFT_COMPLEX STREQUAL "cuFFT")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IComplexBackend.cu
    ComplexProjector.cu
    ComplexIntegrator.cu
  )
set(_use_cuFFT ON)
endif()

if(QUICC_FFT_CHEBYSHEV STREQUAL "cuFFT")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IChebyshevBackend.cu
    GpuDctTools.cu
    ChebyshevProjector.cu
    ChebyshevIntegrator.cu
    ChebyshevEnergy.cu
    ../Fftw/DifferentialSolver.cpp
  )
set(_use_cuFFT ON)
endif()

if(QUICC_FFT_WORLAND STREQUAL "cuFFT")
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    IWorlandBackend.cu
    WorlandProjector.cu
    WorlandIntegrator.cu
  )
set(_use_cuFFT ON)
endif()

if(_use_cuFFT)
quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  CheckCuda.cu
  ICuFftBackend.cu
  GpuMatrix.cu
  GpuCsrMatrix.cu
  GpuPinnedMatrix.cu
  Library.cu
  Tools.cu
)
endif()
