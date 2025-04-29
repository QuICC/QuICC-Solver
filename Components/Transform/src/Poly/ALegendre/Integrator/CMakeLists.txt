quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
  Base/P.cpp
  Base/D1.cpp
  Base/DivLl.cpp
  Base/DivLlD1.cpp
  Base/DivLlDivS1.cpp
  Base/DivLlDivS1Dp.cpp
  Base/DivS1.cpp
  Base/DivS1Dp.cpp
  Base/Ll.cpp
  Base/Ll2.cpp
  Base/LlD1.cpp
  Base/LlDivS1.cpp
  Base/LlDivS1Dp.cpp
  IALegendreIntegrator.cpp
)

if(QUICC_USE_KOKKOS)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    Kokkos/KokkosIALegendreIntegrator.cpp
  )
endif()

if(QUICC_USE_KOKKOS AND NOT QUICC_USE_KOKKOS_CUDA)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    Kokkos/P.cpp
    Kokkos/D1.cpp
    Kokkos/DivLl.cpp
    Kokkos/DivLlD1.cpp
    Kokkos/DivLlDivS1.cpp
    Kokkos/DivLlDivS1Dp.cpp
    Kokkos/DivS1.cpp
    Kokkos/DivS1Dp.cpp
    Kokkos/Ll.cpp
    Kokkos/Ll2.cpp
    Kokkos/LlD1.cpp
    Kokkos/LlDivS1.cpp
    Kokkos/LlDivS1Dp.cpp
  )
endif()

if(QUICC_USE_KOKKOS AND QUICC_USE_KOKKOS_CUDA)
  quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    KokkosCuda/P.cpp
    KokkosCuda/D1.cpp
    KokkosCuda/DivLl.cpp
    KokkosCuda/DivLlD1.cpp
    KokkosCuda/DivLlDivS1.cpp
    KokkosCuda/DivLlDivS1Dp.cpp
    KokkosCuda/DivS1.cpp
    KokkosCuda/DivS1Dp.cpp
    KokkosCuda/Ll.cpp
    KokkosCuda/Ll2.cpp
    KokkosCuda/LlD1.cpp
    KokkosCuda/LlDivS1.cpp
    KokkosCuda/LlDivS1Dp.cpp
  )
endif()
