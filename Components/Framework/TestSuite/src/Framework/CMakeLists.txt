quicc_target_sources(${QUICC_CURRENT_COMPONENT_TEST_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Io.cpp
)

add_subdirectory(Communicators)
add_subdirectory(LoadSplitter)
add_subdirectory(StateFile)
add_subdirectory(Transform)
add_subdirectory(TransformConfigurators)
add_subdirectory(TransformCoordinator)
