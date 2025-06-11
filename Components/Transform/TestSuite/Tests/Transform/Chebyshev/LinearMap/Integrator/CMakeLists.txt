message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)
set(_tests
  P:integrator
  Y1:integrator
  I2:integrator
  I2D1:integrator
  I2D1_I2:integrator
  I2Y1_Zero:integrator
  I2_I2D1:integrator
  I4D1:integrator
  I4D1_I2:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 11
  )
endforeach()

if(QUICC_MULTPRECISION)
  set(_ulp 60)
else()
  set(_ulp 33)
endif()

set(_tests
  I2Y2_Zero:integrator
  I4:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP ${_ulp}
  )
endforeach()

set(_tests
  I2Y1D1Y1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 100
  )
endforeach()

set(_tests
  I2Y2D1Y1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 700
  )
endforeach()

set(_tests
  I4Y3_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 150
  )
endforeach()

set(_tests
  I4Y3D1Y1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    ULP 175
  )
endforeach()
