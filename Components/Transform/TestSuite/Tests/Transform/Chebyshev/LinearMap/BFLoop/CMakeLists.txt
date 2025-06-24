message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)
quicc_clear_test_list()

set(_tests
  P:<base_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
  )
endforeach()

# Merged tests
quicc_add_merged_test(
  ${TestExe}_bfloop_merged
  COMMAND ${TestExe}
  JN 1
)
