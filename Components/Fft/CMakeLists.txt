# CMake setup for Fft backends
set(QUICC_CURRENT_COMPONENT_LIB Fft)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )


add_library(${QUICC_CURRENT_COMPONENT_LIB} INTERFACE)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set(PHS
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Fftw
# option(QUICC_USE_FFTW "Enable Fftw backend" ON)
if(TARGET FFTW::FFTW)
  message(VERBOSE "Fftw backend available")
  add_subdirectory(Fft/Fftw)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE
    QuICC::Fft::Fftw
  )
else()
message(VERBOSE "Fftw backend not available")
endif()

# CuFft
# option(QUICC_USE_CUFFT "Enable CuFft backend" OFF)
find_package(CUDAToolkit)
if(TARGET CUDA::cufft)
  message(VERBOSE "Cufft backend available")
  add_subdirectory(Fft/CuFft)
  target_compile_definitions("${QUICC_CURRENT_COMPONENT_LIB}"
    INTERFACE QUICC_USE_CUFFT)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE
      QuICC::Fft::CuFft
  )
else()
message(VERBOSE "Cufft backend not available")
endif()

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# # Teststuite
# option(QUICC_TESTSUITE_FFT "Enable Fft component testsuite?" OFF)
# if(QUICC_TESTSUITE_FFT)
#   add_subdirectory(TestSuite)
# endif(QUICC_TESTSUITE_FFT)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
