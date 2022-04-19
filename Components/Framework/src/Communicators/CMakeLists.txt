quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Communicator.cpp
  CommunicatorStorage.cpp
  )

add_subdirectory(Converters)
