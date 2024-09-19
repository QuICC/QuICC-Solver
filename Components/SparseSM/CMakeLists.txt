# CMake setup for SparseSM component
set(QUICC_CURRENT_COMPONENT_DIR Components/SparseSM)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_SPARSESM_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    QuICC::Common
  )

# Boost for Bessel
include(QuICCFindBoostRequired)
quicc_find_boost_required()
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
    Boost::headers
)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT SparseSM
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
