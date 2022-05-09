# CMake setup for Framework component
set(QUICC_CURRENT_COMPONENT_DIR Components/Framework)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_FRAMEWORK_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_TRANSFORM_COMPONENT_LIB}
  ${QUICC_POLYNOMIAL_COMPONENT_LIB}
  ${QUICC_PYQUICC_COMPONENT_LIB}
  ${QUICC_IO_COMPONENT_LIB}_ascii
  ${QUICC_IO_COMPONENT_LIB}_config
  ${QUICC_IO_COMPONENT_LIB}_control
  ${QUICC_IO_COMPONENT_LIB}_hdf5
  ${QUICC_IO_COMPONENT_LIB}_xml
  ${QUICC_COMMON_COMPONENT_LIB}
  Profiling::Profiling
  )

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)
