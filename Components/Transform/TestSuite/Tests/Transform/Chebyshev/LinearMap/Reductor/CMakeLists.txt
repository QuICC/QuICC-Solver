message(VERBOSE "Reductors:")

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Energy:reductor
  EnergyY2:reductor
  EnergyD1:reductor
  EnergyD1Y1:reductor
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
  )
endforeach()
