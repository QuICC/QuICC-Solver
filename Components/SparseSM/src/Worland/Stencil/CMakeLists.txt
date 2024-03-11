quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Value.cpp
  ValueDiags.cpp
  D1.cpp
  D1Diags.cpp
  R1D1DivR1.cpp
  R1D1DivR1Diags.cpp
  InsulatingSphere.cpp
  InsulatingSphereDiags.cpp
  ValueD1.cpp
  ValueD1Diags.cpp
  ValueD2.cpp
  ValueD2Diags.cpp
)

add_subdirectory(Chebyshev)
add_subdirectory(SphEnergy)
