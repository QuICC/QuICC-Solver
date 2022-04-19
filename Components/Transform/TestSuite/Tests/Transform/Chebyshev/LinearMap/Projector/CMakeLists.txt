message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)

set(_tests
  P:projector
  DivY1:projector
  DivY2:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 40
  )
endforeach()

set(_tests
  D1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 60
  )
endforeach()

set(_tests
  D4:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 150
  )
endforeach()

set(_tests
  D1Y1:projector
  DivY1D1Y1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 250
  )
endforeach()

set(_tests
  D3:projector
  D2:projector
  SphRadLapl:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 400
  )
endforeach()
