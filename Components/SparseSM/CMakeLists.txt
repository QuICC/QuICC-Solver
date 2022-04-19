# CMake setup for SparseSM component
set(QUICC_CURRENT_COMPONENT_DIR Components/SparseSM)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_SPARSESM_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)
