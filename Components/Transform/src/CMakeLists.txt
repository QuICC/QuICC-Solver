quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ITransformOperator.cpp
  TransformSetup.cpp
)

add_subdirectory(Fft)
add_subdirectory(Poly)
