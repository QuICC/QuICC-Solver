quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ILinearMapOperator.cpp
  DipolarS1.cpp
  QuadrupolarS2.cpp
  R2DivR1F.cpp
  R4DivR1F.cpp
  R4DivR1D1F.cpp
  R2DivR2D1R1F.cpp
  R4DivR2D1R1F.cpp
  R2DivR2FD1R1.cpp
  R4DivR2FD1R1.cpp
)

set(PHS
  DenseSM/Chebyshev/LinearMap/ILineareMapOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
