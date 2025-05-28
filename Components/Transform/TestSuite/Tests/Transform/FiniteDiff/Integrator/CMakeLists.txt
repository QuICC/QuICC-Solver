message(VERBOSE "Integrators:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  FiniteDiff::P:<viewCpu_t>:integrator
  FiniteDiff::P_Zero:<viewCpu_t>:integrator
  FiniteDiff::R1_Zero:<viewCpu_t>:integrator
  FiniteDiff::DivR1_Zero:<viewCpu_t>:integrator
  FiniteDiff::DivR1D1R1_Zero:<viewCpu_t>:integrator
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
    IDS 0 1 2 3 4
  )
endforeach()
