quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IConverter.cpp
  IIndexConv.cpp
  PMIndexConv.cpp
  SHm2lIndexConv.cpp
  SerialConverter.cpp
  SerialConverterBase.cpp
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
