if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    P.cpp
    D1.cpp
    D1R1.cpp
    D1_P.cpp
    DivR1.cpp
    DivR1_Zero.cpp
    DivR1D1R1.cpp
    DivR1D1R1_Zero.cpp
    SphLapl.cpp
    CylLaplh.cpp
    CylLaplh_DivR1D1R1.cpp
    D1CylLaplh.cpp
    D1CylLaplh_D1DivR1D1R1.cpp
    DivR1CylLaplh_Zero.cpp
  )
endif()
