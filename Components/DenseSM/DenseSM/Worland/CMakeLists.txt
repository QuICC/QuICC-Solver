quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IWorlandOperator.cpp
  DivR1F.cpp
  DivR1D1DivR1D1R1F.cpp
  DivR1D1DivR1FD1R1.cpp
  DivR1D1F.cpp
  DivR2D1R1F.cpp
  DivR2FD1R1.cpp
  DivR3D1R1FD1R1.cpp
  DipolarS1.cpp
  QuadrupolarS2.cpp
)

set(PHS
  DenseSM/Worland/IWorlandOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
