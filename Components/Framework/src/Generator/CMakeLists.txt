quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  StateGenerator.cpp
  VisualizationGenerator.cpp
  )

add_subdirectory(States)
add_subdirectory(Visualizers)
