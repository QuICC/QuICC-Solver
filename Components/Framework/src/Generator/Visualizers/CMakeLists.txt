quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ScalarFieldTrivialVisualizer.cpp
  ScalarFieldVisualizer.cpp
  SphericalRadialCylindricalFieldVisualizer.cpp
  SphericalVerticalFieldVisualizer.cpp
  VectorFieldTrivialVisualizer.cpp
  VectorFieldVisualizer.cpp
  )

add_subdirectory(Kernels)
