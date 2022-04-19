# Configure Worland Python interface
configure_file(
  "include/QuICC/Polynomial/Quadrature/WorlandRule.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/Polynomial/Quadrature/WorlandRule.hpp"
  )

# Configure Worland Python interface
configure_file(
  "Python/quicc/geometry/worland/setup.py.in"
  "${PROJECT_BINARY_DIR}/Python/quicc/geometry/worland/setup.py"
  )
