quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IComplexOperator.cpp
  Setup.cpp
  Transform.cpp
)

add_subdirectory(Integrator)
add_subdirectory(Projector)
