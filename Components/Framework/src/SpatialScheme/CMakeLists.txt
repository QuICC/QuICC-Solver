quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Coordinator.cpp
  IBuilder.cpp
  ICosts.cpp
  ICreator.cpp
  ISpatialScheme.cpp
  )

add_subdirectory(1D)
add_subdirectory(2D)
add_subdirectory(3D)
add_subdirectory(Tools)
