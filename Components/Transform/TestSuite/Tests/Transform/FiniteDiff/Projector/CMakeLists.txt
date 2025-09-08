message(VERBOSE "Projectors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create  tests targets
include(ConfigureTests)

set(_tests
  FiniteDiff::P:<viewCpu_t>:projector
  FiniteDiff::R1:<viewCpu_t>:projector
  FiniteDiff::DivR1:<viewCpu_t>:projector
  FiniteDiff::DivR1D1R1:<viewCpu_t>:projector
  FiniteDiff::R1_Zero:<viewCpu_t>:projector
  FiniteDiff::DivR1_Zero:<viewCpu_t>:projector
  FiniteDiff::DivR1D1R1_Zero:<viewCpu_t>:projector
)

quicc_set_ulp(_ulp
  ULP 25 25 25 25
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

set(_tests
  FiniteDiff::D1:<viewCpu_t>:projector
  FiniteDiff::D2:<viewCpu_t>:projector
  FiniteDiff::D3:<viewCpu_t>:projector
  FiniteDiff::D4:<viewCpu_t>:projector
  FiniteDiff::SphLapl:<viewCpu_t>:projector
)

quicc_set_ulp(_ulp
  ULP 600000000
  MPULP 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0
  )
endforeach()

set(_tests
  FiniteDiff::D1:<viewCpu_t>:projector
  FiniteDiff::D2:<viewCpu_t>:projector
  FiniteDiff::D3:<viewCpu_t>:projector
  FiniteDiff::D4:<viewCpu_t>:projector
  FiniteDiff::SphLapl:<viewCpu_t>:projector
)

quicc_set_ulp(_ulp
  ULP 60 60 150
  MPULP 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1 2 3
  )
endforeach()
