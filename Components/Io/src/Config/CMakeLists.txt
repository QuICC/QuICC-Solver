include(ConfigureLibrary)
quicc_add_library(config
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    ConfigurationReader.cpp
    ConfigurationWriter.cpp
    IConfigurationBlock.cpp
    IConfigurationNode.cpp
  CLASSICHEADER
  )

add_subdirectory(Framework)
add_subdirectory(Model)
add_subdirectory(Setup)
add_subdirectory(Simulation)
