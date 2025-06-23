message(VERBOSE "BFLoop:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)
quicc_clear_test_list()

set(_tests
  Poly::P:<base_t>:bfloop
  Poly::P:<viewCpu_t>:bfloop
)
quicc_set_ulp(_ulp
  ULP 160 220 16000 88 130 560 1963 5211 19780
  MPULP 4 4 6 2 2 4 4 8 12
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3 104 105 106 107 108
    ULPS ${_ulp}
  )
endforeach()

# triangular truncations
set(_tests
  Poly::P:<base_t>:bfloop
  # Poly::P:<viewCpu_t>:bfloop
)
quicc_set_ulp(_ulp
  ULP 20 50 148 80 130 560 1963 5211 19780
  MPULP 4 4 6 2 2 4 4 8 12
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 1000 1001 1002 1003 1104 1105 1106 1107
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
    ULPS 5 5 12 5 5 9 11 16 38 76 155
  )
endforeach()

# Merged tests
set(_jN 6)
math(EXPR _jMax "${_jN} - 1")
foreach(_jid RANGE ${_jMax})
  add_test(
    NAME ${TestExe}_bfloop_merged_${_jid}_${_jN}
    COMMAND ${TestExe} --options_file ${CMAKE_CURRENT_BINARY_DIR}/all_tests.txt --jid ${_jid} --jN ${_jN}
    WORKING_DIRECTORY ${QUICC_WORK_DIR}
  )
endforeach()
