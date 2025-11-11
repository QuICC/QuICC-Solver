message(VERBOSE "Reductors:")

# Create reductor tests targets
include(ConfigureTests)
quicc_clear_test_list()

set(_tests
  Energy:<base_t>:reductor
  EnergyY2:<base_t>:reductor
  EnergyD1:<base_t>:reductor
  EnergyD1Y1:<base_t>:reductor
  Energy:<viewCpu_t>:reductor
  EnergyY2:<viewCpu_t>:reductor
  EnergyD1:<viewCpu_t>:reductor
  EnergyD1Y1:<viewCpu_t>:reductor
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
  )
endforeach()

# Merged tests
quicc_add_merged_test(
  ${TestExe}_reductor_merged
  COMMAND ${TestExe}
  JN 4
)
