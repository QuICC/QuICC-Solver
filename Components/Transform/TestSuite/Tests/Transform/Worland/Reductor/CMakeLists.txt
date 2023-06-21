message(VERBOSE "Reductors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Poly::EnergyR2:reductor
  Poly::Energy:reductor
  Poly::EnergySLaplR2:reductor
  Fft::EnergyR2:reductor
  Fft::Energy:reductor
  Fft::EnergySLaplR2:reductor
)
quicc_set_ulp(_ulp
  ULP 30 36 30
  MPULP 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::RadialPower:reductor
  Poly::RadialPowerDivR1:reductor
  Poly::RadialPowerDivR1D1R1:reductor
)
quicc_set_ulp(_ulp
  ULP 60 120 120
  MPULP 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::EnergyD1R1:reductor
  Fft::EnergyD1R1:reductor
)
quicc_set_ulp(_ulp
  ULP 30 201 230
  MPULP 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerD1R1:reductor
  Fft::PowerD1R1:reductor
)
quicc_set_ulp(_ulp
  ULP 30 500 800
  MPULP 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerR2:reductor
  Fft::PowerR2:reductor
)
quicc_set_ulp(_ulp
  ULP 30 300 1600
  MPULP 11 15 70
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::Power:reductor
  Fft::Power:reductor
)
quicc_set_ulp(_ulp
  ULP 30 350 500
  MPULP 11 30 50
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerSLaplR2:reductor
  Fft::PowerSLaplR2:reductor
)
quicc_set_ulp(_ulp
  ULP 30 100 250
  MPULP 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()
