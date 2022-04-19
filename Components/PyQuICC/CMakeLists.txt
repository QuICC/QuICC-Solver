# CMake setup for PyQuICC component
set(QUICC_CURRENT_COMPONENT_DIR Components/PyQuICC)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_PYQUICC_COMPONENT_LIB})

# Set path to installed python package and scripts
set(QUICC_PYTHON_DIR "${PROJECT_BINARY_DIR}/Python")

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Python code
include(cmake.d/setup/Python.cmake)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  include/
  ${PROJECT_BINARY_DIR}/include
  )

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC 
  ${QUICC_COMMON_COMPONENT_LIB}
  Python::Python Python::NumPy
  )

# Install python files
add_dependencies(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CURRENT_COMPONENT_LIB}_updatepy)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)
