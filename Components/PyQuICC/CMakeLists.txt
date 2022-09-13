# CMake setup for PyQuICC component
set(QUICC_CURRENT_COMPONENT_DIR Components/PyQuICC)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_PYQUICC_COMPONENT_LIB})

# Set path to installed python package and scripts
set(QUICC_PYTHON_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python")

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Python code
include(cmake.d/setup/Python.cmake)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC 
  QuICC::Common
  Python::Python
  Python::NumPy
  )

# Install python files
add_dependencies(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CURRENT_COMPONENT_LIB}_updatepy)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT PyQuICC
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
