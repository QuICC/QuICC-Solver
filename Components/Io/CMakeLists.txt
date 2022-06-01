# CMake setup for Io component
set(QUICC_CURRENT_COMPONENT_DIR Components/Io)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_IO_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library for IO
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  hdf5::hdf5
)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Io
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
