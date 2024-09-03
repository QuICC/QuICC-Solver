message(VERBOSE "Projectors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create  tests targets
include(ConfigureTests)

set(_tests
  Value::P:<base_t>:projector
  Insulating::P:<base_t>:projector
  NoSlip::P:<base_t>:projector
)

quicc_set_ulp(_ulp
  ULP 8 94 55 60
  MPULP 1 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Value::D1:<base_t>:projector
  Value::DivR1_Zero:<base_t>:projector
  Value::DivR1D1R1_Zero:<base_t>:projector
  Value::SphLapl:<base_t>:projector
  Insulating::D1:<base_t>:projector
  Insulating::DivR1_Zero:<base_t>:projector
  Insulating::DivR1D1R1_Zero:<base_t>:projector
  Insulating::SphLapl:<base_t>:projector
  NoSlip::D1:<base_t>:projector
  NoSlip::DivR1_Zero:<base_t>:projector
  NoSlip::DivR1D1R1_Zero:<base_t>:projector
  NoSlip::SphLapl:<base_t>:projector
)

quicc_set_ulp(_ulp
  ULP 9 60 49 51
  MPULP 1 1 2 2
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

# Split projector tests

set(_split_tests
  Value::P:<base_t>:projector
  Insulating::P:<base_t>:projector
  NoSlip::P:<base_t>:projector
)

foreach(_test ${_split_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS  104 104 108 108 108 108 108 108
    ULPS 150 150 15000 15000 15000 15000 15000 15000
    # node:rank, node*12:rank, node*36:rank
    SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
  )
endforeach()

# Perf tests

set(_perf_tests
  Value::P:<base_t>:projector
  Insulating::P:<base_t>:projector
  NoSlip::P:<base_t>:projector
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
    list(APPEND _perf_ulps 15000)
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
