message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)
set(_tests
  Poly::P:<base_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 20
  )
endforeach()
