# Configure the Python embedding wrapper
configure_file(
  "${QUICC_PYQUICC_DIR}/include/QuICC/PyQuICC/Config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/PyQuICC/Config.hpp"
  )

# Update python files for PyQuICC
add_custom_target(${QUICC_CURRENT_COMPONENT_LIB}_updatepy)
add_custom_command(TARGET ${QUICC_CURRENT_COMPONENT_LIB}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  "${CMAKE_CURRENT_SOURCE_DIR}/Python"
  "${CMAKE_BINARY_DIR}/Python"
  COMMENT "Copying Python files for PyQuICC"
  VERBATIM
  )
