# CMake setup for DenseOp
set(QUICC_CURRENT_COMPONENT_LIB DenseOp)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

add_library(${QUICC_CURRENT_COMPONENT_LIB} INTERFACE)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

add_subdirectory(DenseOp/Worland)

set(PHS
  DenseOp/IDenseOpOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Dependencies
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} INTERFACE
  QuICC::Polynomial
  ${QUICC_SPARSESM_COMPONENT_LIB}
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# # Teststuite
# option(QUICC_TESTSUITE_Cuda "Enable Cuda component testsuite?" OFF)
# if(QUICC_TESTSUITE_Cuda)
#   add_subdirectory(TestSuite)
# endif(QUICC_TESTSUITE_Cuda)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
)
