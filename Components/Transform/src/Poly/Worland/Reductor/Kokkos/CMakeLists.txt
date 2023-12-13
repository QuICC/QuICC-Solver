if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB_KK} ${QUICC_CMAKE_SRC_VISIBILITY}
    KokkosIWorlandReductor.cpp
    KokkosIWorlandPower.cpp
    KokkosIWorlandRadialPower.cpp
  )
endif()

