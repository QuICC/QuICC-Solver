quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ILinearMapOperator.cpp
  RpDivR1F.cpp
  R0F.cpp
  RpFD1.cpp
  RpFD2.cpp
  RpFD3.cpp
  RpFD4.cpp
  RpDivR1FC.cpp
  RpDivR1CF.cpp
  RpDivR1D1F.cpp
  RpDivR1D1FC.cpp
  RpDivR1D1CF.cpp
  RpDivR2D1R1F.cpp
  RpDivR2D1R1FC.cpp
  RpDivR2FD1R1.cpp
  RpDivR2CFD1R1.cpp
  RpDivR1D1DivR1FD1R1.cpp
  RpDivR1D1DivR1D1R1F.cpp
  RpDivR3D1R1FD1R1.cpp
  IqRpDivR1F.cpp
  IqRpDivR1FC.cpp
  IqRpDivR1CF.cpp
  IqRpDivR1D1F.cpp
  IqRpDivR1D1FC.cpp
  IqRpDivR1D1CF.cpp
  IqRpDivR2D1R1F.cpp
  IqRpDivR2D1R1FC.cpp
  IqRpDivR2FD1R1.cpp
  IqRpDivR2CFD1R1.cpp
  IqRpDivR1D1DivR1FD1R1.cpp
  IqRpDivR1D1DivR1D1R1F.cpp
  IqRpDivR3D1R1FD1R1.cpp
  IProjCrossOperator.cpp
  ProjCurlTorCrossTor.cpp
  ProjCurlTorCrossPol.cpp
  ProjCurlPolCrossTor.cpp
  ProjCurlPolCrossPol.cpp
  ProjCurlCPolCrossTor.cpp
  ProjCurlCPolCrossPol.cpp
  ProjCurlCurlTorCrossTor.cpp
  ProjCurlCurlTorCrossPol.cpp
  ProjCurlCurlPolCrossTor.cpp
  ProjCurlCurlPolCrossPol.cpp
  ProjCurlCurlCPolCrossTor.cpp
  ProjCurlCurlCPolCrossPol.cpp
  RadialTorPolFunction.cpp
)

set(PHS
  DenseSM/Chebyshev/LinearMap/ProjoCurlTorCrossTor.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlTorCrossPol.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlPolCrossTor.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlPolCrossPol.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlCurlTorCrossTor.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlCurlTorCrossPol.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlCurlPolCrossTor.hpp
  DenseSM/Chebyshev/LinearMap/ProjoCurlCurlPolCrossPol.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

add_subdirectory(Utils)
