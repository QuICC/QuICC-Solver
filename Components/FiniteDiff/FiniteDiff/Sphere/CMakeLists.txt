target_sources(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
  Operator.cpp
  Id.cpp
  Integral.cpp
  IntegralR2.cpp
  D1.cpp
  D1R1.cpp
  Overr1.cpp
  Overr1D1R1.cpp
  R1.cpp
  SLapl.cpp
  UniformRadialGrid.cpp
  )

add_subdirectory(Boundary)
