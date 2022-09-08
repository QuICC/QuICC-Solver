include(ConfigureLibrary)
quicc_add_library(ascii
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    AsciiFile.cpp
    IAsciiReader.cpp
    IAsciiWriter.cpp
    DirectAsciiWriter.cpp
    StdOutPipe.cpp
    CflTags.cpp
    CflWriter.cpp
  CLASSICHEADER
)
