include(ConfigureLibrary)
quicc_add_library(hdf5
  PARENT ${QUICC_CURRENT_COMPONENT_LIB}
  LINK PUBLIC QuICC::Common hdf5::hdf5
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  SRC
    Hdf5File.cpp
    IHdf5Reader.cpp
    IHdf5Writer.cpp
    IHdf5NWriter.cpp
  CLASSICHEADER
)
