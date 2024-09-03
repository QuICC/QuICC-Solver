message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Value::P:<base_t>:bfloop
  Insulating::P:<base_t>:bfloop
  NoSlip::P:<base_t>:bfloop
)
quicc_set_ulp(_ulp
  ULP 160 220 16000 88 130 560 1963 5211 19780
  MPULP 4 4 6 2 2 4 4 8 12
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS 0 1 2 3 104 105 106 107 108
    ULPS ${_ulp}
  )
endforeach()
