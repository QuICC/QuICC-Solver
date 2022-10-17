message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Poly::P:bfloop
  )
quicc_set_ulp(_ulp
  ULP 16 13 63 39 130 478 1963 3870 7625
  MPULP 3 2 2 2 2 2 3 6 12
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2 3 204 205 206 207 208
    ULPS ${_ulp}
  )
endforeach()

# Create BFLoop tests targets
set(_tests
  Fft::P:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0 1 2 3 204 205 206 207 208 209 210
    ULPS 2 3 5 3 4 6 8 16 21 28 41
  )
endforeach()
