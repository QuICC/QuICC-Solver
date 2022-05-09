# CMake setup for Transform component
set(QUICC_CURRENT_COMPONENT_DIR Components/Transform)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_TRANSFORM_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Components/Transform/cmake.d")

# Setup for the different transforms
include(setup/ALegendre)
include(setup/Worland)
include(setup/SparseLinearAlgebra)
include(setup/FFT)
include(setup/CUDA)


# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  ${QUICC_SPARSESM_COMPONENT_LIB}
  ${QUICC_POLYNOMIAL_COMPONENT_LIB}
  Profiling::Profiling
  BLAS::BLAS # banded/tridiag solve/mul TODO: remove dep
  )

if(QUICC_WORLAND_BACKEND STREQUAL "BLAS")
  if(NOT TARGET BLAS::BLAS)
    message(SEND_ERROR "${QUICC_CURRENT_COMPONENT_LIB} requires BLAS, but it was not found.")
  endif()
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
    BLAS::BLAS
    )
endif()

include(MatchAny)
match_any(NAME "QUICC_FFT_" STRING "FFTW")
if(FFTW_IS_USED)
  message(DEBUG "FFTW is used")
  if(NOT FFTW_FOUND)
    message(SEND_ERROR "${QUICC_CURRENT_COMPONENT_LIB} requires FFTW, but it was not found.")
  endif()
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE FFTW::FFTW)
  quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE OPTION QUICC_FFTPLAN)
endif()
# match_any cuFFT


# Add source directory
add_subdirectory(src)

# Teststuite
include(setup/TestSuite)
