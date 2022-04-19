message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)
set(_tests
  Mixed::P:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
  )
endforeach()
