# CMake setup for DenseSM
set(QUICC_CURRENT_COMPONENT_DIR Components/DenseSM)
set(QUICC_CURRENT_COMPONENT_LIB DenseSM)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  QuICC::Common
  QuICC::Polynomial
  ${QUICC_SPARSESM_COMPONENT_LIB}
  ${QUICC_PYQUICC_COMPONENT_LIB}
)

# Add source directory
add_subdirectory(DenseSM)

set(PHS
  DenseSM/IDenseSMOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
# option(QUICC_TESTSUITE_Cuda "Enable Cuda component testsuite?" OFF)
include(cmake.d/setup/TestSuite.cmake)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
)
