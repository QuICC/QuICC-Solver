# CMake setup for View component
set(QUICC_CURRENT_COMPONENT_LIB ViewOps)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

# components
add_subdirectory(ViewOps/Blas)
add_subdirectory(ViewOps/Pointwise)
add_subdirectory(ViewOps/Reduction)
add_subdirectory(ViewOps/Fourier)
add_subdirectory(ViewOps/Quadrature)
add_subdirectory(ViewOps/ALegendre)
add_subdirectory(ViewOps/Worland)

add_library(${QUICC_CURRENT_COMPONENT_LIB} INTERFACE)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Add public headers
set(PHS
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# View
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE
  QuICC::View
  QuICC::Operator
  QuICC::ViewOps::Fourier
  QuICC::ViewOps::ALegendre
  QuICC::ViewOps::Worland
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
option(QUICC_TESTSUITE_VIEWOPS "Enable ViewOps component testsuite?" OFF)
if(QUICC_TESTSUITE_VIEWOPS)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_VIEWOPS)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
)
