quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  JacobiBase.cpp
)

if(NOT QUICC_MULTPRECISION)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
    JacobiAsymptotic.cpp
)
endif()
