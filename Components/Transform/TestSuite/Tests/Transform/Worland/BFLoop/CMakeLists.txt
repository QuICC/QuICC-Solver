message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Poly::P:<base_t>:bfloop
  )
quicc_set_ulp(_ulp
  ULP 160 220 16000 88 20 50 130 80 130 560 1963 5211 19780 130 560 1963 5211 19780
  MPULP 4 4 6 2 4 4 6 2 2 4 4 8 12 2 4 4 8 12
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3 1000 1001 1002 1003 104 105 106 107 108 1104 1105 1106 1107
    1108
    ULPS ${_ulp}
  )
endforeach()

# Create BFLoop tests targets
set(_tests
  Fft::P:<base_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3 104 105 106 107 108 109 110
    ULPS 5 5 12 5 5 6 11 16 21 28 41
  )
endforeach()
