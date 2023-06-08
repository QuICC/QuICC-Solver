message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)
set(_tests
  Mixed::P:<base_t>:bfloop
  Mixed::P:<viewCpu_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3 4
    ULPS 11 11 11 11 11
  )
endforeach()
