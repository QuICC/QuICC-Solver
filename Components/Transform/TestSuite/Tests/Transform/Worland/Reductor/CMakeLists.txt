message(VERBOSE "Reductors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergyD1R1:reductor
  Poly::EnergyR2:reductor
  Poly::Energy:reductor
  Poly::EnergySLAPLR2:reductor
  )
quicc_set_ulp(_ulp ULP 30 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergyR2:reductor
  Poly::Energy:reductor
  Poly::EnergySLAPLR2:reductor
  )
quicc_set_ulp(_ulp ULP 30 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergyD1R1:reductor
  )
quicc_set_ulp(_ulp ULP 190 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergyD1R1:reductor
  )
quicc_set_ulp(_ulp ULP 220 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergySLAPLR2:reductor
  Poly::EnergyR2:reductor
  Poly::Energy:reductor
  )
quicc_set_ulp(_ulp ULP 30 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerD1R1:reductor
  Poly::PowerR2:reductor
  Poly::Power:reductor
  Poly::PowerSLAPLR2:reductor
  )
quicc_set_ulp(_ulp ULP 30 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerD1R1:reductor
  )
quicc_set_ulp(_ulp ULP 500 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerR2:reductor
  )
quicc_set_ulp(_ulp ULP 300 MPULP 15)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::Power:reductor
  )
quicc_set_ulp(_ulp ULP 350 MPULP 30)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerSLAPLR2:reductor
  )
quicc_set_ulp(_ulp ULP 100 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerD1R1:reductor
  )
quicc_set_ulp(_ulp ULP 800 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerR2:reductor
  )
quicc_set_ulp(_ulp ULP 1600 MPULP 70)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::Power:reductor
  )
quicc_set_ulp(_ulp ULP 500 MPULP 50)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::PowerSLAPLR2:reductor
  )
quicc_set_ulp(_ulp ULP 250 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()
