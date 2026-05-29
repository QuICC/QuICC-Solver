quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IWorlandOperator.cpp
  DivR1F.cpp
  DivR1CF.cpp
  DivR1FC.cpp
  DivR1D1DivR1D1R1F.cpp
  DivR1D1DivR1FD1R1.cpp
  DivR1D1FD1R1.cpp
  DivR1D1FD1R1C.cpp
  DivR1D1F.cpp
  DivR1D1CF.cpp
  DivR1D1FC.cpp
  DivR2F.cpp
  DivR2D1R1F.cpp
  DivR1FD1.cpp
  DivR1FD1R1.cpp
  DivR2FD1R1.cpp
  DivR2FD3R1.cpp
  DivR2CFD1R1.cpp
  DivR2D1R1FC.cpp
  DivR3FD2R1.cpp
  DivR3D1R1FD1R1.cpp
  F.cpp
  FC.cpp
  FC2.cpp
  IProjCrossOperator.cpp
  ProjCurlCPolCrossPol.cpp
  ProjCurlCPolCrossTor.cpp
  ProjCurlCurlCPolCrossPol.cpp
  ProjCurlCurlCPolCrossTor.cpp
  ProjCurlCurlPolCrossPol.cpp
  ProjCurlCurlPolCrossTor.cpp
  ProjCurlCurlTorCrossPol.cpp
  ProjCurlCurlTorCrossTor.cpp
  ProjCurlPolCrossPol.cpp
  ProjCurlPolCrossTor.cpp
  ProjCurlTorCrossPol.cpp
  ProjCurlTorCrossTor.cpp
  ProjPolViscNu.cpp
  ProjPolViscD1.cpp
  ProjPolViscD2.cpp
  ProjPolViscD3.cpp
  ProjPolViscD4.cpp
  ProjTorViscNu.cpp
  ProjTorViscD1.cpp
  ProjTorViscD2.cpp
  Tools.cpp
)

set(PHS
  DenseSM/Worland/IWorlandOperator.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
