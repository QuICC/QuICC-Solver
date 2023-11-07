if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    Energy.cpp
    EnergyR2.cpp
    EnergyD1R1.cpp
    EnergySLaplR2.cpp
    Power.cpp
    PowerR2.cpp
    PowerD1R1.cpp
    PowerSLaplR2.cpp
    RadialPower.cpp
    RadialPowerDivR1.cpp
    RadialPowerDivR1D1R1.cpp
  )
endif()

