message(STATUS "... adding simulation tools")

# CMake setup for Simulation tools
set(QUICC_CURRENT_TOOL Simulation)
set(QUICC_CURRENT_TOOL_DIR Tools/${QUICC_CURRENT_TOOL})
set(QUICC_CURRENT_TOOL_LIB Tools_${QUICC_CURRENT_TOOL})
set(QUICC_CURRENT_TOOL_LIB_NAMESPACE "Tools::${QUICC_CURRENT_TOOL}::")

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${QUICC_CURRENT_TOOL_DIR}/cmake.d")

# Teststuite
include(setup/ModelTestSuite)

add_library(${QUICC_CURRENT_TOOL_LIB} INTERFACE)

target_include_directories(${QUICC_CURRENT_TOOL_LIB}
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Add library alias (as the imported one)
add_library(${QUICC_NAMESPACE}Tools::${QUICC_CURRENT_TOOL}
  ALIAS "${QUICC_CURRENT_TOOL_LIB}")

# Export info
quicc_export_target(${QUICC_CURRENT_TOOL_LIB}
  COMPONENT ${QUICC_CURRENT_TOOL_LIB}
)
