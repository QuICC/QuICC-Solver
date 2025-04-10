message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)
include(SetUlp)

quicc_set_ulp(_ulps
  ULP 11 11 11 11 11 11 20 25
  MPULP 11 11 11 11 11 11 20 25
  )

set(_tests
  P:<base_t>:bfloop
  P:<viewCpu_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 3 4 5 6 7
  )
endforeach()
