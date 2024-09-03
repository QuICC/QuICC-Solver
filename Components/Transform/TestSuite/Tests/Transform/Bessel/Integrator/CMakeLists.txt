message(VERBOSE "Integrators:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Value::P:<base_t>:integrator
  Value::P_Zero:<base_t>:integrator
  Value::R1_Zero:<base_t>:integrator
  Insulating::P:<base_t>:integrator
  Insulating::P_Zero:<base_t>:integrator
  Insulating::R1_Zero:<base_t>:integrator
  NoSlip::P:<base_t>:integrator
  NoSlip::P_Zero:<base_t>:integrator
  NoSlip::R1_Zero:<base_t>:integrator
)

quicc_set_ulp(_ulp
  ULP 16 7 115 80 78
  MPULP 3 2 2 2 2
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()


set(_tests
  Value::DivR1_Zero:<base_t>:integrator
  Insulating::DivR1_Zero:<base_t>:integrator
  NoSlip::DivR1_Zero:<base_t>:integrator
)

quicc_set_ulp(_ulp
  ULP 32 51 226 172 191 32 42 148 75 73
  MPULP 2 2 2 2 2 3 3 3 3 3
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Value::DivR1D1R1_Zero:<base_t>:integrator
  Insulating::DivR1D1R1_Zero:<base_t>:integrator
  NoSlip::DivR1D1R1_Zero:<base_t>:integrator
)

quicc_set_ulp(_ulp
  ULP 453 94 118 119 122 605 97 119 141 134
  MPULP 6 2 1 1 1 9 2 1 2 1
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()


#-----------------
# Database tests
#-----------------

set(_tests_templated
  Value::P:<base_t>:integrator
  Insulating::P:<base_t>:integrator
  NoSlip::P:<base_t>:integrator
)

quicc_set_ulp(_ulp
  ULP 111 346 1331 2240 12000
  MPULP 8 20 12 18 20
  )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 104 105 106 107 108
  )
endforeach()

# Split integrator tests

set(_split_tests
  Value::P:<base_t>:integrator
  Insulating::P:<base_t>:integrator
  NoSlip::P:<base_t>:integrator
)

foreach(_test ${_split_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS  104 104 108 108 108 108 108 108
    ULPS 120 120 13000 11000 11000 12000 11000 11000
    # node:rank, node*12:rank, node*36:rank
    SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
  )
endforeach()

# Perf tests

set(_perf_tests
  Value::P:<base_t>:integrator
  Insulating::P:<base_t>:integrator
  NoSlip::P:<base_t>:integrator
)

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
      KEYWORD Boundary::Operator:Template:type
      IDS  ${_perf_ids}
      ULPS ${_perf_ulps}
      SPLITS ${_perf_splits}
      STEPS ${_perf_steps}
      PERFONLY
    )
  endforeach()
endforeach()
