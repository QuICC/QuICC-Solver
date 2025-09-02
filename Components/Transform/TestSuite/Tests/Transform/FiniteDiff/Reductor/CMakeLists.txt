message(VERBOSE "Reductors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  FiniteDiff::Energy:<viewCpu_t>:reductor
  FiniteDiff::EnergyR2:<viewCpu_t>:reductor
  FiniteDiff::EnergySLaplR2:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 36 30 32
  MPULP 11 11 11 11
)
set(_ids
  0 1 2 3
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  FiniteDiff::EnergyD1R1:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 201 230 240
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()
