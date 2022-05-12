# Configure the Python embedding wrapper
configure_file(
  "${QUICC_PYQUICC_DIR}/include/QuICC/PyQuICC/Config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/PyQuICC/Config.hpp"
  )

# Update python files for PyQuICC
# akin to install, but run everytime a py file is changed
add_custom_target(${QUICC_CURRENT_COMPONENT_LIB}_updatepy)
add_custom_command(TARGET ${QUICC_CURRENT_COMPONENT_LIB}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  "${CMAKE_CURRENT_SOURCE_DIR}/Python"
  "${QUICC_PYTHON_DIR}"
  COMMENT "Copying Python files for PyQuICC"
  VERBATIM
  )
