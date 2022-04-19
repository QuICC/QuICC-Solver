message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)
set(_tests
  Poly::DivLlDivS1DpTest:integrator
  Poly::DivLlDivS1Test:integrator
  Poly::DivLl:integrator
  Poly::DivS1Dp:integrator
  Poly::DivS1:integrator
  Poly::LlDivS1Dp:integrator
  Poly::LlDivS1:integrator
  Poly::Ll:integrator
  Poly::P:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 100
  )
endforeach()

set(_tests
  Poly::D1Test:integrator
  Poly::LlD1:integrator
  Poly::Ll:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 150
  )
endforeach()

set(_tests
  Poly::Ll2:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1400
  )
endforeach()

# Special integrator tests
set(_tests
  Poly::P:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 104 105 106 107 108
    ULP 100
  )
endforeach()

# Special integrator tests
set(_tests
  Poly::P:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 109
    ULP 15000
  )
endforeach()

# Special integrator tests
set(_tests
  Poly::P:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 110
    ULP 80000
  )
endforeach()
