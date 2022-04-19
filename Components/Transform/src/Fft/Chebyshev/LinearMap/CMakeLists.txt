quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Transform.cpp
)

add_subdirectory(Projector)
add_subdirectory(Integrator)
add_subdirectory(Reductor)
