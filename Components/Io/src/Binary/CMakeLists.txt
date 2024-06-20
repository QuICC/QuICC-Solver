include(ConfigureLibrary)
quicc_add_library(binary
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common QuICC::Environment
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    BinaryFile.cpp
    IBinaryWriter.cpp
    DirectBinaryWriter.cpp
  CLASSICHEADER
)
