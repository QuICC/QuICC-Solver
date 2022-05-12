# CMake setup for Polynomial component
set(QUICC_CURRENT_COMPONENT_DIR Components/Polynomial)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_POLYNOMIAL_COMPONENT_LIB})

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
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

add_subdirectory(Bessel)
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  Polynomial::Bessel
  )

# Additional configuration
include(setup/Worland_configure)

# Add source directory
add_subdirectory(src)

# Testsuite
include(setup/TestSuite)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Polynomial
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
