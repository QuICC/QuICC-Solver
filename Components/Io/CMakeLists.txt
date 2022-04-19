# CMake setup for Io component
set(QUICC_CURRENT_COMPONENT_DIR Components/Io)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_IO_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library for ASCII
add_library(${QUICC_CURRENT_COMPONENT_LIB}_ascii "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_ascii PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_ascii PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_ascii PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Create component library for binary
add_library(${QUICC_CURRENT_COMPONENT_LIB}_binary "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_binary PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_binary PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_binary PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Create component library for config
add_library(${QUICC_CURRENT_COMPONENT_LIB}_config "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_config PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_config PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_config PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Create component library for control
add_library(${QUICC_CURRENT_COMPONENT_LIB}_control "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_control PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_control PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_control PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Create component library for HDF5
add_library(${QUICC_CURRENT_COMPONENT_LIB}_hdf5 "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_hdf5 PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_hdf5 PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_hdf5 PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  hdf5::hdf5
  )

# Create component library for XML
add_library(${QUICC_CURRENT_COMPONENT_LIB}_xml "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_xml PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_xml PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_xml PUBLIC
  ${QUICC_COMMON_COMPONENT_LIB}
  )

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)
