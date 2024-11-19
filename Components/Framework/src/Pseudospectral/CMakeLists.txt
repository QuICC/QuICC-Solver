quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Coordinator.cpp
  AddGraph.cpp
  ComputeNonLinear.cpp
  Utils.cpp
)

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    QuICC::Graph
)
