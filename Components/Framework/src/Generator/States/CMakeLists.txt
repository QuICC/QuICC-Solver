quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  RandomScalarState.cpp
  RandomStateData.cpp
  RandomVectorState.cpp
  CartesianExactScalarState.cpp
  CartesianExactVectorState.cpp
  ShellExactScalarState.cpp
  ShellExactVectorState.cpp
  SphereExactScalarState.cpp
  SphereExactVectorState.cpp
  )

add_subdirectory(Kernels)
