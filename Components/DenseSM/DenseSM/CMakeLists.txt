quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IEmbeddedSMOperator.cpp
  IMatrixSMOperator.cpp
)

add_subdirectory(Bessel)
add_subdirectory(Worland)
