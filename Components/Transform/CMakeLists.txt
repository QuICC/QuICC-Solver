# CMake setup for Transform component
set(QUICC_CURRENT_COMPONENT_DIR Components/Transform)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_TRANSFORM_COMPONENT_LIB}_nonkokkos)
set(QUICC_CURRENT_COMPONENT_LIB_KK ${QUICC_TRANSFORM_COMPONENT_LIB}_kokkos)

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Components/Transform/cmake.d")

# Setup for the different transforms
include(setup/Worland)
include(setup/SparseLinearAlgebra)
include(setup/FFT)
include(setup/CUDA)


# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )

include(setup/TruncateQI)

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  QuICC::Types
  QuICC::Environment
  ${QUICC_SPARSESM_COMPONENT_LIB}
  QuICC::Polynomial
  QuICC::Profiling
  QuICC::ViewOps
  QuICC::Memory
  BLAS::BLAS # banded/tridiag solve/mul TODO: remove dep
  )

# Wrappers depends on kk part
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
    ${QUICC_CURRENT_COMPONENT_LIB_KK}
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

if(QUICC_USE_KOKKOS)
  include(setup/KokkosAddLibrary)
  quicc_add_library_kokkos(${QUICC_CURRENT_COMPONENT_LIB_KK} "")
else()
  add_library(${QUICC_CURRENT_COMPONENT_LIB_KK} "")
endif()

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB_KK}
  PUBLIC
    ${QUICC_CURRENT_COMPONENT_LIB}
    QuICC::KokkosWrapper
)

# VkFFT
if(TARGET VkFFT AND QUICC_USE_VKFFT)
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE VkFFT)
endif()

# PfSolve
if(TARGET PfSolve AND QUICC_USE_PFSOLVE)
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE PfSolve)
endif()

# View Operator wrappers
add_subdirectory(Wrappers)

# Add source directory
add_subdirectory(src)

if(QUICC_USE_KOKKOS)
  # fix file properties
  quicc_fix_target_sources_kokkos(${QUICC_CURRENT_COMPONENT_LIB_KK})
endif()

# Group components
add_library(${QUICC_TRANSFORM_COMPONENT_LIB} INTERFACE)
target_link_libraries(${QUICC_TRANSFORM_COMPONENT_LIB} INTERFACE
  ${QUICC_CURRENT_COMPONENT_LIB}
  ${QUICC_CURRENT_COMPONENT_LIB_KK})
add_library(QuICC::Transform ALIAS ${QUICC_TRANSFORM_COMPONENT_LIB})

# Teststuite
include(setup/TestSuite)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB_KK}
  COMPONENT TransformKokkos
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )

quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Transform
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )

quicc_export_target(${QUICC_TRANSFORM_COMPONENT_LIB}
  COMPONENT TransformAll
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )

