quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ILinearMapOperator.cpp
  DipolarS1.cpp
  QuadrupolarS2.cpp
  R2DivR1F.cpp
  R4DivR1F.cpp
  R4DivR1FC.cpp
  R4DivR1CF.cpp
  R4DivR1D1F.cpp
  R4DivR1D1FC.cpp
  R4DivR1D1CF.cpp
  R2DivR2D1R1F.cpp
  R4DivR2D1R1F.cpp
  R4DivR2D1R1FC.cpp
  R2DivR2FD1R1.cpp
  R4DivR2FD1R1.cpp
  R4DivR2CFD1R1.cpp
  R4DivR1D1DivR1FD1R1.cpp
  R4DivR1D1DivR1D1R1F.cpp
  R4DivR3D1R1FD1R1.cpp
  #IProjCrossOperator.cpp
  #ProjCurlTorCrossTor.cpp
  #ProjCurlTorCrossPol.cpp
  #ProjCurlPolCrossTor.cpp
  #ProjCurlPolCrossPol.cpp
  #ProjCurlCPolCrossTor.cpp
  #ProjCurlCPolCrossPol.cpp
  #ProjCurlCurlTorCrossTor.cpp
  #ProjCurlCurlTorCrossPol.cpp
  #ProjCurlCurlPolCrossTor.cpp
  #ProjCurlCurlPolCrossPol.cpp
  #ProjCurlCurlCPolCrossTor.cpp
  #ProjCurlCurlCPolCrossPol.cpp
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
