quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
  D1.cpp
  DivLl.cpp
  DivLlD1.cpp
  DivLlDivS1.cpp
  DivLlDivS1Dp.cpp
  DivS1.cpp
  DivS1Dp.cpp
  IALegendreIntegrator.cpp
  Ll2.cpp
  Ll.cpp
  LlD1.cpp
  LlDivS1.cpp
  LlDivS1Dp.cpp
  P.cpp
)

if(QUICC_USE_KOKKOS)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    PIALegendreIntegrator.cpp
    PP.cpp
    Cuda/PP.cpp
    Cuda/PD1.cpp
    Cuda/PLl.cpp
    Cuda/PLlDivS1.cpp
    Cuda/PLlDivS1Dp.cpp
    Cuda/PLl2.cpp
    Cuda/PDivLl.cpp
    Cuda/PLlD1.cpp
    Cuda/PDivLlD1.cpp
    Cuda/PDivS1.cpp
    Cuda/PDivS1Dp.cpp
    Cuda/PDivLlDivS1.cpp
    Cuda/PDivLlDivS1Dp.cpp
  )
endif(QUICC_USE_KOKKOS)
