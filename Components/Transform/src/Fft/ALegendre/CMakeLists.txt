quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IALegendreOperator.cpp
  Setup.cpp
  Tools.cpp
  Transform.cpp
)

add_subdirectory(Projector)
add_subdirectory(Integrator)
