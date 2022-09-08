include(ConfigureLibrary)
quicc_add_library(control
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    ControlFile.cpp
    ControlInterface.cpp
  CLASSICHEADER
)
