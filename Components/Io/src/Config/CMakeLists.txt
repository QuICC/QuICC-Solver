quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB}_config ${QUICC_CMAKE_SRC_VISIBILITY}
  ConfigurationReader.cpp
  ConfigurationWriter.cpp
  IConfigurationBlock.cpp
  IConfigurationNode.cpp
  )

add_subdirectory(Framework)
add_subdirectory(Model)
add_subdirectory(Setup)
add_subdirectory(Simulation)
