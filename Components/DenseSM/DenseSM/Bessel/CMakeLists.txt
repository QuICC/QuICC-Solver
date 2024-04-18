quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  CoriolisQm.cpp
  CoriolisQp.cpp
  Geostrophic2Tor.cpp
  ProjFitEnergyR2.cpp
  ProjFitEnergyPol.cpp
  Tor2Geostrophic.cpp
  Tor2GridS.cpp
)

add_subdirectory(details)
