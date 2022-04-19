message(VERBOSE "Reductors:")

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Poly::EnergyD1R1:reductor
  Poly::EnergyR2:reductor
  Poly::Energy:reductor
  Poly::EnergySLAPLR2:reductor
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
  )
endforeach()
