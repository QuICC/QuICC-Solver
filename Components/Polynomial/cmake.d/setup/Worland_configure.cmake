# Configure Worland header
configure_file(
  "include/QuICC/Polynomial/Quadrature/WorlandRule.hpp.in"
  "${CMAKE_BINARY_DIR}/include/QuICC/Polynomial/Quadrature/WorlandRule.hpp"
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
      PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Polynomial/Quadrature/WorlandRule.hpp"
)

# Configure Worland Python interface
configure_file(
  "Python/quicc/geometry/worland/setup.py.in"
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Python/quicc/geometry/worland/setup.py"
)

# Update python files for PyQuICC
# akin to install, but run everytime a py file is changed
add_custom_target(${QUICC_CURRENT_COMPONENT_LIB}_updatepy)
add_custom_command(TARGET ${QUICC_CURRENT_COMPONENT_LIB}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Python/quicc/geometry/worland/setup.py"
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python/quicc/geometry/worland/setup.py"
  COMMENT "Copying Python files for PyQuICC"
  VERBATIM
)
