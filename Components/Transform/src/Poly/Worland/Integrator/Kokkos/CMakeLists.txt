if(QUICC_USE_KOKKOS)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
      KokkosIWorlandIntegrator.cpp
  )
endif()


if(QUICC_USE_KOKKOS AND NOT QUICC_USE_KOKKOS_CUDA)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    P.cpp
    P_Zero.cpp
    R1_Zero.cpp
    DivR1_Zero.cpp
    DivR1D1R1_Zero.cpp
    I2.cpp
    I2_Zero.cpp
    I2DivR1_Zero.cpp
    I2DivR1D1R1_Zero.cpp
    I4DivR1_Zero.cpp
    I4DivR1D1R1_I2.cpp
    I4DivR1D1R1_Zero.cpp
    I6DivR1_Zero.cpp
  )
endif()
