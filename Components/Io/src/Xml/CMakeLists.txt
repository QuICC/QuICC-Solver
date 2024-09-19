include(ConfigureLibrary)
quicc_add_library(xml
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common QuICC::Environment External::rxml
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    XmlFile.cpp
    IXmlReader.cpp
    IXmlWriter.cpp
  CLASSICHEADER
)
