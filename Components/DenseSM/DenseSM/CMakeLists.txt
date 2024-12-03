quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IMatrixSMOperator.cpp
)

set(PHS
  DenseSM/IMatrixSMOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

add_subdirectory(Chebyshev)
add_subdirectory(Worland)
