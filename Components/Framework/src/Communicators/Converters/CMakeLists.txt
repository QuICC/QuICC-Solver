quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IConverter.cpp
  IIndexConv.cpp
  NoIndexConv.cpp
  PMIndexConv.cpp
  SerialConverter.cpp
  SerialConverterBase.cpp
  SHlIndexConv.cpp
  SHmIndexConv.cpp
  )

if(QUICC_MPI)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    MpiConverter.cpp
    MpiConverterAllToAll.cpp
    MpiConverterBase.cpp
    MpiConverterSendRecv.cpp
    MpiConverterTools.cpp
    )
endif(QUICC_MPI)
