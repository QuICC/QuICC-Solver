# CMake setup for Io component
set(QUICC_CURRENT_COMPONENT_DIR Components/Io)
set(QUICC_CURRENT_COMPONENT_LIB Io)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "Io::" )

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Inteface target for all components
add_library("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE)
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}"
  INTERFACE
    QuICC::Io::ascii
    QuICC::Io::binary
    QuICC::Io::config
    QuICC::Io::control
    QuICC::Io::hdf5
    QuICC::Io::xml
    )

# Alias
add_library(${QUICC_NAMESPACE}Io ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Io
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
