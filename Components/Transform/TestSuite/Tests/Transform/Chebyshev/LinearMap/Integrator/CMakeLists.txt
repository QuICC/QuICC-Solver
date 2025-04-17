message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)
set(_tests
  P:<base_t>:integrator
  P_Zero:<base_t>:integrator
  Y1:<base_t>:integrator
  Y1_Zero:<base_t>:integrator
  I2:<base_t>:integrator
  I2D1:<base_t>:integrator
  I2D1_I2:<base_t>:integrator
  I2Y1_Zero:<base_t>:integrator
  I2_I2D1:<base_t>:integrator
  I4D1:<base_t>:integrator
  I4D1_I2:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP 11
    IDS 0 1 2 3 4 5
  )
endforeach()

if(QUICC_MULTPRECISION)
  set(_ulp 60)
else()
  set(_ulp 33)
endif()

set(_tests
  I2Y2_Zero:<base_t>:integrator
  I4:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP ${_ulp}
    IDS 0 1 2 3 4 5
  )
endforeach()

set(_tests
  I2Y1D1Y1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP 100
    IDS 0 1 2 3 4 5
  )
endforeach()

set(_tests
  I2Y2D1Y1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP 700
    IDS 0 1 2 3 4 5
  )
endforeach()

set(_tests
  I4Y3_Zero:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP 150
    IDS 0 1 2 3 4 5
  )
endforeach()

set(_tests
  I4Y3D1Y1_Zero:<base_t>:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP 175
    IDS 0 1 2 3 4 5
  )
endforeach()
