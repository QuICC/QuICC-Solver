quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ITransform.cpp
  ALegendreTransform.cpp
  CartesianChebyshevTransform.cpp
  ComplexFourierTransform.cpp
  CylinderWorlandTransform.cpp
  MixedFourierTransform.cpp
  ShellChebyshevTransform.cpp
  SphereFftWorlandTransform.cpp
  SphereWorlandTransform.cpp
  AnnulusChebyshevTransform.cpp
)

add_subdirectory(Backward)
add_subdirectory(Forward)
add_subdirectory(Reductor)
