target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IWorlandOperator.cpp
  ProjFitEnergyR2.cpp
  ProjFitEnergyPol.cpp
  I2ProjFitEnergyR2.cpp
  I2ProjFitEnergyPol.cpp
  IGeostrophicOperator.cpp
  Geostrophic2Tor.cpp
  GeostrophicEnergy.cpp
  PyGeostrophic2Tor.cpp
  GeostrophicAngularMomentum.cpp
  GeostrophicSolidBody.cpp
  Tools.cpp
  Tor2Geostrophic.cpp
  Tor2EGeostrophic.cpp
  Tor2GridS.cpp
  Tor2Weights.cpp
)

add_subdirectory(details)
