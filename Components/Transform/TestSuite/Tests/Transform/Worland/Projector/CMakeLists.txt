message(VERBOSE "Projectors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create  tests targets
include(ConfigureTests)

set(_tests
  Poly::P:<base_t>:projector
  Poly::P:<viewCpu_t>:projector
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::P:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::P:<kokkos_t>:projector
    )
endif()

quicc_set_ulp(_ulp
  ULP 8 94 55 60
  MPULP 1 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()


# triangular truncation
set(_tests
  Poly::P:<base_t>:projector
  # Poly::P:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 8 140 55 60
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003
  )
endforeach()


# \todo DivR1 test missing

set(_tests
  Poly::D1:<base_t>:projector
  Poly::D1R1:<base_t>:projector
  Poly::DivR1_Zero:<base_t>:projector
  Poly::DivR1D1R1:<base_t>:projector
  Poly::DivR1D1R1_Zero:<base_t>:projector
  Poly::SphLapl:<base_t>:projector
  Poly::D1:<viewCpu_t>:projector
  Poly::D1R1:<viewCpu_t>:projector
  Poly::DivR1_Zero:<viewCpu_t>:projector
  Poly::DivR1D1R1:<viewCpu_t>:projector
  Poly::DivR1D1R1_Zero:<viewCpu_t>:projector
  Poly::SphLapl:<viewCpu_t>:projector
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::D1:<viewGpu_t>:projector
    Poly::D1R1:<viewGpu_t>:projector
    Poly::DivR1_Zero:<viewGpu_t>:projector
    Poly::DivR1D1R1:<viewGpu_t>:projector
    Poly::DivR1D1R1_Zero:<viewGpu_t>:projector
    Poly::SphLapl:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_KOKKOS)
list(APPEND _tests
  Poly::D1:<kokkos_t>:projector
  Poly::D1R1:<kokkos_t>:projector
  Poly::DivR1_Zero:<kokkos_t>:projector
  Poly::DivR1D1R1:<kokkos_t>:projector
  Poly::DivR1D1R1_Zero:<kokkos_t>:projector
  Poly::SphLapl:<kokkos_t>:projector
)
endif()

quicc_set_ulp(_ulp
  ULP 9 60 49 51
  MPULP 1 1 2 2
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

# triangular truncation
set(_tests
  Poly::D1:<base_t>:projector
  Poly::D1R1:<base_t>:projector
  Poly::DivR1_Zero:<base_t>:projector
  Poly::DivR1D1R1:<base_t>:projector
  Poly::DivR1D1R1_Zero:<base_t>:projector
  Poly::SphLapl:<base_t>:projector
  # Poly::D1:<viewCpu_t>:projector
  # Poly::D1R1:<viewCpu_t>:projector
  # Poly::DivR1_Zero:<viewCpu_t>:projector
  # Poly::DivR1D1R1:<viewCpu_t>:projector
  # Poly::DivR1D1R1_Zero:<viewCpu_t>:projector
  # Poly::SphLapl:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 9 60 49 51
  MPULP 1 1 2 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003
  )
endforeach()

set(_tests
  Poly::CylLaplh:<base_t>:projector
  Poly::CylLaplh_DivR1D1R1:<base_t>:projector
  Poly::D1_P:<base_t>:projector
  Poly::D1CylLaplh:<base_t>:projector
  Poly::D1CylLaplh_D1DivR1D1R1:<base_t>:projector
  Poly::DivR1:<base_t>:projector
  Poly::DivR1CylLaplh_Zero:<base_t>:projector
  Poly::CylLaplh:<viewCpu_t>:projector
  Poly::CylLaplh_DivR1D1R1:<viewCpu_t>:projector
  Poly::D1_P:<viewCpu_t>:projector
  Poly::D1CylLaplh:<viewCpu_t>:projector
  Poly::D1CylLaplh_D1DivR1D1R1:<viewCpu_t>:projector
  Poly::DivR1:<viewCpu_t>:projector
  Poly::DivR1CylLaplh_Zero:<viewCpu_t>:projector
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::CylLaplh:<viewGpu_t>:projector
    Poly::CylLaplh_DivR1D1R1:<viewGpu_t>:projector
    Poly::D1_P:<viewGpu_t>:projector
    Poly::D1CylLaplh:<viewGpu_t>:projector
    Poly::D1CylLaplh_D1DivR1D1R1:<viewGpu_t>:projector
    Poly::DivR1:<viewGpu_t>:projector
    Poly::DivR1CylLaplh_Zero:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_KOKKOS)
    set(_tests
      Poly::CylLaplh:<kokkos_t>:projector
      Poly::CylLaplh_DivR1D1R1:<kokkos_t>:projector
      Poly::D1CylLaplh:<kokkos_t>:projector
      Poly::D1_P:<kokkos_t>:projector
      Poly::D1CylLaplh_D1DivR1D1R1:<kokkos_t>:projector
      Poly::DivR1CylLaplh_Zero:<kokkos_t>:projector
    )
quicc_set_ulp(_ulp
  ULP 1 1 1 1
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
    DISABLED ON
  )
endforeach()
endif()

# triangular truncation
set(_tests
  Poly::CylLaplh:<base_t>:projector
  Poly::CylLaplh_DivR1D1R1:<base_t>:projector
  Poly::D1_P:<base_t>:projector
  Poly::D1CylLaplh:<base_t>:projector
  Poly::D1CylLaplh_D1DivR1D1R1:<base_t>:projector
  Poly::DivR1:<base_t>:projector
  Poly::DivR1CylLaplh_Zero:<base_t>:projector
  # Poly::CylLaplh:<viewCpu_t>:projector
  # Poly::CylLaplh_DivR1D1R1:<viewCpu_t>:projector
  # Poly::D1_P:<viewCpu_t>:projector
  # Poly::D1CylLaplh:<viewCpu_t>:projector
  # Poly::D1CylLaplh_D1DivR1D1R1:<viewCpu_t>:projector
  # Poly::DivR1:<viewCpu_t>:projector
  # Poly::DivR1CylLaplh_Zero:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 1 1 1 1
  MPULP  1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003
    DISABLED ON
  )
endforeach()


######################
# FFT based algorithm
######################

set(_tests
  Fft::P:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 17 5 5
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1_Zero:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 6 5 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::D1:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 4 4 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1D1R1_Zero:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 4 3 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::SphLapl:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 1 4 5 3
    IDS  0 1 2 3
  )
endforeach()

# Split projector tests

set(_split_tests
  Fft::P:<base_t>:projector
  Poly::P:<base_t>:projector
  Poly::P:<viewCpu_t>:projector
)
if(TARGET QuICC::Cuda)
    list(APPEND _split_tests
        Poly::P:<viewGpu_t>:projector
    )
endif()
if(QUICC_USE_KOKKOS)
  list(APPEND _split_tests
      Poly::P:<kokkos_t>:projector
  )
endif()

foreach(_test ${_split_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  104 104 108 108 108 108 108 108
    ULPS 150 150 15000 15000 15000 15000 15000 15000
    # node:rank, node*12:rank, node*36:rank
    SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
  )
endforeach()

# Perf tests

set(_perf_tests
  Poly::P:<base_t>:projector
  Poly::P:<viewCpu_t>:projector
  Fft::P:<base_t>:projector
)
if(TARGET QuICC::Cuda)
    list(APPEND _perf_tests
        Poly::P:<viewGpu_t>:projector
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _perf_tests
        Poly::P:<kokkos_t>:projector
        )
endif()

foreach(_perf_id 108 109)
  if(${_perf_id} EQUAL 109)
    set(_perf_step_slow 50)
    set(_perf_step_fast 50)
  else()
    set(_perf_step_slow 500)
    set(_perf_step_fast 1000)
  endif()

  set(_perf_splits_slow 8:0 16:0 32:0 64:0 96:0 128:0 192:0 256:0 288:0 384:0)
  set(_perf_splits_fast 576:0 768:0 1152:0 1536:0 2304:0 3072:0 4608:0 9216:0)

  set(_perf_ids )
  set(_perf_ulps )
  set(_perf_steps )
  set(_perf_splits )
  foreach(_s ${_perf_splits_slow})
    list(APPEND _perf_ids ${_perf_id})
    list(APPEND _perf_ulps 15000)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_slow})
  endforeach()
  foreach(_s ${_perf_splits_fast})
    list(APPEND _perf_ids ${_perf_id})
    list(APPEND _perf_ulps 15000)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_fast})
  endforeach()

  foreach(_test ${_perf_tests})
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Algorithm::Operator:Template:type
      IDS  ${_perf_ids}
      ULPS ${_perf_ulps}
      SPLITS ${_perf_splits}
      STEPS ${_perf_steps}
      PERFONLY
      )
  endforeach()
endforeach()
