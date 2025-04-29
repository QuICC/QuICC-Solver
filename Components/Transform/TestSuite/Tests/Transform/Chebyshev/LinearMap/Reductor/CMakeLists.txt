message(VERBOSE "Reductors:")

# Create reductor tests targets
include(ConfigureTests)
set(_tests
  Energy:<base_t>:reductor
  EnergyY2:<base_t>:reductor
  EnergyD1:<base_t>:reductor
  EnergyD1Y1:<base_t>:reductor
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
  )
endforeach()
