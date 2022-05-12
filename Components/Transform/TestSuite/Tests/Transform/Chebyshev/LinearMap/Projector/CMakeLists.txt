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

if(QUICC_MULTPRECISION)
  set(_ulp 110)
else()
  set(_ulp 60)
endif()

set(_tests
  D1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP ${_ulp}
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

if(QUICC_MULTPRECISION)
  set(_ulp 900)
else()
  set(_ulp 250)
endif()

set(_tests
  D1Y1:projector
  DivY1D1Y1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP ${_ulp}
  )
endforeach()

if(QUICC_MULTPRECISION)
  set(_ulp 600)
else()
  set(_ulp 400)
endif()

set(_tests
  D3:projector
  D2:projector
  SphRadLapl:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP ${_ulp}
  )
endforeach()
