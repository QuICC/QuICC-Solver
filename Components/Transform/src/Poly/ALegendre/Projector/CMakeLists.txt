quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
  D1.cpp
  DivS1.cpp
  DivS1D1S1.cpp
  DivS1Dp.cpp
  IALegendreProjector.cpp
  Ll.cpp
  LlD1.cpp
  LlDivS1.cpp
  LlDivS1Dp.cpp
  P.cpp
)


if(QUICC_USE_KOKKOS)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    PIALegendreProjector.cpp
    PP.cpp
    Cuda/PP.cpp
    Cuda/PD1.cpp
    Cuda/PDivS1.cpp
    Cuda/PDivS1Dp.cpp
    Cuda/PLl.cpp
    Cuda/PLlD1.cpp
    Cuda/PLlDivS1.cpp
    Cuda/PLlDivS1Dp.cpp
  )
endif(QUICC_USE_KOKKOS)
