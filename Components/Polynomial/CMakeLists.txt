# CMake setup for Polynomial component
set(QUICC_CURRENT_COMPONENT_DIR Components/Polynomial)
set(QUICC_CURRENT_COMPONENT_LIB Polynomial)

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Components/Polynomial/cmake.d")

# Setup for the different polynomials
include(setup/ALegendre)
include(setup/Jacobi)
include(setup/Legendre)
include(setup/Worland)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)

set_property(TARGET ${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

# Components
add_subdirectory(Bessel)
add_subdirectory(SphericalBessel)
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  QuICC::Types
  Polynomial::Bessel
  Polynomial::SphericalBessel
)
# PfSolve
if(TARGET PfSolve AND QUICC_USE_PFSOLVE)
  target_compile_definitions("${QUICC_CURRENT_COMPONENT_LIB}"
    INTERFACE QUICC_USE_PFSOLVE)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE
    PfSolve
  )
  target_compile_definitions("${QUICC_CURRENT_COMPONENT_LIB}"
      INTERFACE "QUICC_HAS_CUDA_BACKEND")
else()
  message(VERBOSE "PfSolve backend not available")
endif()
# Additional configuration
include(setup/Worland_configure)

# Install python files
add_dependencies(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CURRENT_COMPONENT_LIB}_updatepy)
if(QUICC_CURRENT_UPDATEPY_TARGET)
   add_dependencies(${QUICC_CURRENT_UPDATEPY_TARGET} ${QUICC_CURRENT_COMPONENT_LIB}_updatepy)
endif()
set(QUICC_CURRENT_UPDATEPY_TARGET "${QUICC_CURRENT_COMPONENT_LIB}_updatepy" CACHE STRING "Make dependencies across updatepy" FORCE)

# Add source directory
add_subdirectory(src)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# Testsuite
include(setup/TestSuite)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
