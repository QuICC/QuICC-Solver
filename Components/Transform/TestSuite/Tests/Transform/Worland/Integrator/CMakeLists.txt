message(VERBOSE "Integrators:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

#set(_tests
#  Poly::I4DivR1D1R1_I2:integrator
#  Poly::I6CylLaplh_I4D1R1:integrator
#  Poly::I6DivR1D1R1_I4:integrator
#  )

set(_tests
  Poly::P:<base_t>:integrator
  Poly::P_Zero:<base_t>:integrator
  Poly::R1_Zero:<base_t>:integrator
  Poly::P:<viewCpu_t>:integrator
  Poly::P_Zero:<viewCpu_t>:integrator
  Poly::R1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::P:<viewGpu_t>:integrator
        Poly::P_Zero:<viewGpu_t>:integrator
        Poly::R1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::P:<kokkos_t>:integrator
        Poly::P_Zero:<kokkos_t>:integrator
        Poly::R1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 16 7 115 80 78
  MPULP 3 2 2 2 2
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

# triangular truncation tests
set(_tests
  Poly::P:<base_t>:integrator
  Poly::P_Zero:<base_t>:integrator
  Poly::R1_Zero:<base_t>:integrator
  # Poly::P:<viewCpu_t>:integrator
  # Poly::P_Zero:<viewCpu_t>:integrator
  # Poly::R1_Zero:<viewCpu_t>:integrator
)
quicc_set_ulp(_ulp
  ULP 16 10 115 80 78
  MPULP 3 2 2 2 2
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003 1004
  )
endforeach()


set(_tests
  Poly::DivR1_Zero:<base_t>:integrator
  Poly::DivR1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::DivR1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::DivR1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 32 51 226 172 191 32 42 148 75 73
  MPULP 2 2 2 2 2 3 3 3 3 3
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I2:<base_t>:integrator
  Poly::I2_Zero:<base_t>:integrator
  Poly::I2:<viewCpu_t>:integrator
  Poly::I2_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I2:<viewGpu_t>:integrator
        Poly::I2_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I2:<kokkos_t>:integrator
        Poly::I2_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 14 20 170 59 36
  MPULP 1 4 8 4 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

# triangular truncation tests
set(_tests
  Poly::I2:<base_t>:integrator
  Poly::I2_Zero:<base_t>:integrator
  # Poly::I2:<viewCpu_t>:integrator
  # Poly::I2_Zero:<viewCpu_t>:integrator
)
quicc_set_ulp(_ulp
  ULP 18 20 170 59 120
  MPULP 2 8 26 10 5
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003 1004
  )
endforeach()


set(_tests
  Poly::I2DivR1_Zero:<base_t>:integrator
  Poly::I2DivR1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I2DivR1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I2DivR1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 40 47 145 73 86 41 206 191 41 51
  MPULP 2 3 14 4 3 1 10 11 2 3
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I4DivR1_Zero:<base_t>:integrator
  Poly::I4DivR1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I4DivR1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I4DivR1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 77 127 666 55 91 74 615 453 48 70
  MPULP 2 11 121 9 4 2 127 93 4 6
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I6DivR1_Zero:<base_t>:integrator
  Poly::I6DivR1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I6DivR1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I6DivR1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 53 160 1578 76 51 48 846 1258 66 104
  MPULP 2 11 147 6 8 2 165 92 2 10
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::DivR1D1R1_Zero:<base_t>:integrator
  Poly::DivR1D1R1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::DivR1D1R1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::DivR1D1R1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 453 94 118 119 122 605 97 119 141 134
  MPULP 6 2 1 1 1 9 2 1 2 1
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I2DivR1D1R1_Zero:<base_t>:integrator
  Poly::I2DivR1D1R1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I2DivR1D1R1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I2DivR1D1R1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 30 383 1340 1104 1196 32 214 1017 374 925
  MPULP 2 17 64 42 40 2 10 44 31 24
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I4DivR1D1R1_Zero:<base_t>:integrator
  Poly::I4DivR1D1R1_Zero:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests
        Poly::I4DivR1D1R1_Zero:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests
        Poly::I4DivR1D1R1_Zero:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 73 1262 10831 2375 2005 73 699 6655 1747 2191
  MPULP 2 79 464 121 55 2 58 328 49 51
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests_templated
  Poly::P:<base_t>:integrator
  Poly::P:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _tests_templated
        Poly::P:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _tests_templated
        Poly::P:<kokkos_t>:integrator
    )
endif()

quicc_set_ulp(_ulp
  ULP 111 346 1331 2240 12000
  MPULP 8 20 12 18 20
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 104 105 106 107 108
  )
endforeach()

# triangular truncation
set(_tests_templated
  Poly::P:<base_t>:integrator
  # Poly::P:<viewCpu_t>:integrator
)
quicc_set_ulp(_ulp
  ULP 111 346 1331 3500 12000
  MPULP 2 2 3 4 6
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1104 1105 1106 1107 1108
  )
endforeach()


set(_tests_templated
  Poly::I2:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 68 150 476 1176 180 180 476 1176
  MPULP 8 6 8 9 8 6 8 9
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 104 105 106 107 1104 1105 1106 1107
  )
endforeach()

######################
# FFT based algorithm
######################

set(_tests_templated
  Fft::P:<base_t>:integrator
  Fft::P_Zero:<base_t>:integrator
  Fft::R1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 4 4 6 5 5
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests_templated
  Fft::DivR1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 5 20 11 11 5 5 8 6 5
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests_templated
  Fft::I2:<base_t>:integrator
  Fft::I2_Zero:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 18 22 176 65 31
  MPULP 6 6 55 19 9
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests_templated
  Fft::I2DivR1_Zero:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 16 23 189 48 78 13 213 205 40 59
  MPULP 10 5 33 13 19 7 54 76 9 7
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests_templated
  Fft::I4DivR1_Zero:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 55 23 732 59 51 46 602 582 44 22
  MPULP 17 16 102 9 19 11 209 280 11 11
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests_templated
  Fft::DivR1D1R1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 16 6 5 4 4 18 8 3 4 4
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests_templated
  Fft::I2DivR1D1R1_Zero:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 11 353 1359 1067 1258 11 193 1111 393 884
  MPULP 4 75 278 95 149 3 57 107 163 97
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests_templated
  Fft::I4DivR1D1R1_Zero:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 44 1213 10190 2217 1803 38 658 6820 1591 2111
  MPULP 6 142 550 400 236 4 175 520 58 180
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests_templated
  Fft::P:<base_t>:integrator
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 8 16 14 21 25
    IDS 104 105 106 107 108
  )
endforeach()

set(_tests_templated
  Fft::I2:<base_t>:integrator
  )
quicc_set_ulp(_ulp
  ULP 56 130 472 1145
  MPULP 74 20 26 33
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 104 105 106 107
  )
endforeach()

# Split integrator tests

set(_split_tests
  Fft::P:<base_t>:integrator
  Poly::P:<base_t>:integrator
  Poly::P:<viewCpu_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _split_tests
        Poly::P:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _split_tests
        Poly::P:<kokkos_t>:integrator
        )
endif()

foreach(_test ${_split_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  104 104 108 108 108 108 108 108
    ULPS 120 120 13000 11000 11000 12000 11000 11000
    # node:rank, node*12:rank, node*36:rank
    SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
  )
endforeach()

# Perf tests

set(_perf_tests
  Poly::P:<base_t>:integrator
  Poly::P:<viewCpu_t>:integrator
  Fft::P:<base_t>:integrator
)
if(TARGET QuICC::Cuda)
    list(APPEND _perf_tests
        Poly::P:<viewGpu_t>:integrator
    )
endif()
if(QUICC_USE_KOKKOS)
    list(APPEND _perf_tests
        Poly::P:<kokkos_t>:integrator
        )
endif()

foreach(_perf_id 108 109)
  if(${_perf_id} EQUAL 109)
    set(_perf_step_slow 10)
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
    list(APPEND _perf_ulps 20800)
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
