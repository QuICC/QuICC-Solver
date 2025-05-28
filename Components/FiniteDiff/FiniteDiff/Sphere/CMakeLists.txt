target_sources(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
  Operator.cpp
  Id.cpp
  D1.cpp
  Overr1.cpp
  Overr1D1R1.cpp
  R1.cpp
  SLapl.cpp
  UniformRadialGrid.cpp
  )

add_subdirectory(Boundary)
