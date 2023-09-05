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

# Add option to force independent MPIO for regular data
set( _regularIndep "QUICC_MPIO_REGULAR_FORCE_INDEPENDENT")
option(${_regularIndep} "Force independent MPIO for regular data" OFF)
mark_as_advanced(${_regularIndep})
if(${_regularIndep})
  target_compile_definitions(Io_hdf5 PUBLIC ${_regularIndep})
endif()

# Add option to force independent MPIO for irregular data
set( _irregularIndep "QUICC_MPIO_IRREGULAR_FORCE_INDEPENDENT")
option(${_irregularIndep} "Force independent MPIO for irregular data" OFF)
mark_as_advanced(${_irregularIndep})
if(${_irregularIndep})
  target_compile_definitions(Io_hdf5 PUBLIC ${_irregularIndep})
endif()
