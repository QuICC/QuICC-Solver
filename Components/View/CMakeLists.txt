# CMake setup for View component
set(QUICC_CURRENT_COMPONENT_LIB View)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

add_library(${QUICC_CURRENT_COMPONENT_LIB} INTERFACE)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set(PHS
  "View/View.hpp"
  "View/ViewBase.hpp"
  "View/ViewDense.hpp"
  "View/ViewSparse.hpp"
  "View/Attributes.hpp"
  "View/ViewUtils.hpp"
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Span patch
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE
  QuICC::Std
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
option(QUICC_TESTSUITE_VIEW "Enable View component testsuite?" OFF)
if(QUICC_TESTSUITE_VIEW)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_VIEW)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
    )
