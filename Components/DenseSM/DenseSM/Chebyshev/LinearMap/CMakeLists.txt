quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ILinearMapOperator.cpp
  DipolarS1.cpp
  QuadrupolarS2.cpp
  DivR1F.cpp
)

set(PHS
  DenseSM/Chebyshev/LinearMap/ILineareMapOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
