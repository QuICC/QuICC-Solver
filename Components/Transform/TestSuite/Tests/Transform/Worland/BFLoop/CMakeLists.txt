message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create BFLoop tests targets
include(ConfigureTests)

set(_tests
  Poly::P:bfloop
  )
quicc_set_ulp(_ulp ULP 160 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 0
    ULP ${_ulp}
  )
endforeach()

set(_tests
  Poly::P:bfloop
  )
quicc_set_ulp(_ulp ULP 220 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 1
    ULP ${_ulp}
  )
endforeach()

set(_tests
  Poly::P:bfloop
  )
quicc_set_ulp(_ulp ULP 16000 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 2
    ULP ${_ulp}
  )
endforeach()

# Create BFLoop tests targets
include(ConfigureTests)
set(_tests
  Fft::P:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 11
  )
endforeach()
