quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
  Setup.cpp
  Transform.cpp
  IALegendreOperator.cpp
)

add_subdirectory(Integrator)
add_subdirectory(Projector)
