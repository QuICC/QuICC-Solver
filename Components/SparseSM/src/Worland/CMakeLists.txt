quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IDiags.cpp
  I2.cpp
  I2Diags.cpp
  I4.cpp
  I4Diags.cpp
  I4D1R1.cpp
  I4D1R1Diags.cpp
  I6.cpp
  I6Diags.cpp
  I6CylLaplh.cpp
  I6CylLaplhDiags.cpp
  R2.cpp
  R2Diags.cpp
)

add_subdirectory(Chebyshev)
add_subdirectory(Legendre)
add_subdirectory(CylEnergy)
add_subdirectory(SphEnergy)
