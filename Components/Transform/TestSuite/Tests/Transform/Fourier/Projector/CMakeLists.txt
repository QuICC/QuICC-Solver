message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Mixed::P:projector
  Mixed::D1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    ULP 25
  )
endforeach()

set(_tests
  Mixed::D2:projector
  Mixed::D3:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    ULP 60
  )
endforeach()

set(_tests
  Mixed::D4:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    ULP 250
  )
endforeach()

set(_tests
  Complex::P:projector
  Complex::D1:projector
  Complex::D2:projector
  Complex::Lapl2D:projector
  Complex::Ds1Lapl2D:projector
  Complex::Mean:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
  )
endforeach()

set(_tests
  Complex::D4:projector
  Complex::Df1Lapl2D:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    ULP 25
  )
endforeach()

set(_tests
  Complex::D3:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    ULP 50
  )
endforeach()
