message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Poly::P:bfloop
  )
quicc_set_ulp(_ulp
  ULP 16 18 153 88 130 560 1963 5211 19780
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
    ULPS 5 5 6 5 5 6 11 16 21 28 41
  )
endforeach()
