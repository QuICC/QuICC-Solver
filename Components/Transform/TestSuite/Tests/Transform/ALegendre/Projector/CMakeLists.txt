message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::P:projector
  Poly::D1:projector
  Poly::Ll:projector
  Poly::LlDivS1:projector
  Poly::LlDivS1Dp:projector
  Poly::LlD1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 100
  )
endforeach()

set(_tests
  Poly::DivS1:projector
  Poly::DivS1Dp:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 200
  )
endforeach()

# Special projector tests
set(_tests
  Poly::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 104 105 106 107 108
    ULP 100
  )
endforeach()

# Special projector tests
set(_tests
  Poly::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 109
    ULP 150
  )
endforeach()

# Special projector tests
set(_tests
  Poly::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 110
    ULP 500
  )
endforeach()
