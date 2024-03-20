quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IEmbeddedOperator.cpp
  ProjFitEnergyR2.cpp
  ProjFitEnergyPol.cpp
  I2ProjFitEnergyR2.cpp
  I2ProjFitEnergyPol.cpp
  IGeostrophicOperator.cpp
  Geostrophic2Tor.cpp
  PyGeostrophic2Tor.cpp
  GeostrophicAngularMomentum.cpp
  Tools.cpp
  Tor2Geostrophic.cpp
  Tor2EGeostrophic.cpp
  Tor2GridS.cpp
)

set(PHS
  DenseSM/Worland/GetSuperN.hpp
  DenseSM/Worland/Operator.hpp
  DenseSM/Worland/OperatorWithMean.hpp
  DenseSM/Worland/OperatorIN.hpp
  DenseSM/Worland/OperatorINWithMean.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
