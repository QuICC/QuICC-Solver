quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Coordinator.cpp
  ComputeNonLinear.cpp
)

if (QUICC_USE_MLIR_GRAPH)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    Utils.cpp
    AddGraph.cpp
  )
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
    PUBLIC
      QuICC::Graph
  )
endif()