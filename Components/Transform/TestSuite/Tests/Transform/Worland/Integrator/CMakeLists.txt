message(VERBOSE "Integrators:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create integrator tests targets
include(ConfigureTests)
#set(_tests
#  Poly::I4DivR1D1R1_I2Test:integrator
#  Poly::I6CylLaplh_I4D1R1Test:integrator
#  Poly::I6DivR1D1R1_I4Test:integrator
#  )

set(_tests
  Poly::R1Test:integrator
  Poly::PTest:integrator
  Poly::DivR1Test:integrator
  )
quicc_set_ulp(_ulp ULP 100 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 0 1 2
  )
endforeach()

set(_tests
  Poly::DivR1Test:integrator
  )
quicc_set_ulp(_ulp ULP 100 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3
  )
endforeach()

set(_tests
  Poly::R1Test:integrator
  Poly::PTest:integrator
  )
quicc_set_ulp(_ulp ULP 220 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3
  )
endforeach()

set(_tests
  Poly::R1Test:integrator
  Poly::PTest:integrator
  Poly::DivR1Test:integrator
  )
quicc_set_ulp(_ulp ULP 600 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 4
  )
endforeach()

set(_tests
  Poly::I2_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 45 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 0 1
  )
endforeach()

set(_tests
  Poly::I2_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 2500 MPULP 50)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3
  )
endforeach()

set(_tests
  Poly::I2_ZeroTest:integrator
  )
#quicc_set_ulp(_ulp ULP 5800 MPULP 100)
quicc_set_ulp(_ulp ULP 1 MPULP 1)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 4
  )
endforeach()

set(_tests
  Poly::I2_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 20000 MPULP 350)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 2
  )
endforeach()

set(_tests
  Poly::I2DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 40 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 0 1
  )
endforeach()

set(_tests
  Poly::I4DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 40 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 0
  )
endforeach()

set(_tests
  Poly::I4DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 1300 MPULP 60)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 1
  )
endforeach()

set(_tests
  Poly::I6DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 55 MPULP 11)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 0
  )
endforeach()

set(_tests
  Poly::I6DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 4500 MPULP 300)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 1
  )
endforeach()

set(_tests
  Poly::I2DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 1100 MPULP 100)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3 4
  )
endforeach()

set(_tests
  Poly::I2DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 20000 MPULP 500)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 2
  )
endforeach()

set(_tests
  Poly::I4DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 600 MPULP 100)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3
  )
endforeach()

set(_tests
  Poly::I4DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 800 MPULP 100)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 4
  )
endforeach()

set(_tests
  Poly::I4DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 2500 MPULP 70)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 2
  )
endforeach()

set(_tests
  Poly::I6DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 360 MPULP 20)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 3 4
  )
endforeach()

set(_tests
  Poly::I6DivR1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 2500 MPULP 350)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP ${_ulp}
    IDS 2
  )
endforeach()

set(_tests
  Poly::DivR1D1R1Test:integrator
  Poly::I2DivR1D1R1_ZeroTest:integrator
  Poly::I4DivR1D1R1_ZeroTest:integrator
  )
quicc_set_ulp(_ulp ULP 100000 MPULP 2500)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
  )
endforeach()

# Disable Poly tests temporarily
set_tests_properties(
  TransformWorlandTests_Poly_DivR1D1R1_integrator
  TransformWorlandTests_Poly_I2_Zero_integrator_id4_ulp1
  TransformWorlandTests_Poly_I2DivR1D1R1_Zero_integrator
  TransformWorlandTests_Poly_I4DivR1D1R1_Zero_integrator
  PROPERTIES DISABLED TRUE
  )

set(_tests
  Fft::R1Test:integrator
  Fft::PTest:integrator
  Fft::DivR1Test:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 21
    IDS 0 1 2
  )
endforeach()

set(_tests
  Fft::R1Test:integrator
  Fft::PTest:integrator
  Fft::DivR1Test:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 45
    IDS 3
  )
endforeach()

set(_tests
  Fft::R1Test:integrator
  Fft::PTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 11
    IDS 4
  )
endforeach()

set(_tests
  Fft::I2_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 50
    IDS 0 1
  )
endforeach()

set(_tests
  Fft::I2_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1800
    IDS 3
  )
endforeach()

set(_tests
  Fft::I2_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 3500
    IDS 4
  )
endforeach()

set(_tests
  Fft::I2_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 20000
    IDS 2
  )
endforeach()

set(_tests
  Fft::I2DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 50
    IDS 0 1
  )
endforeach()

set(_tests
  Fft::I4DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 350
    IDS 0 1
  )
endforeach()

set(_tests
  Fft::I2DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 25000
    IDS 2
  )
endforeach()

set(_tests
  Fft::I2DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1100
    IDS 3
  )
endforeach()

set(_tests
  Fft::I2DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 600
    IDS 4
  )
endforeach()

set(_tests
  Fft::I4DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 550
    IDS 3
  )
endforeach()

set(_tests
  Fft::I4DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1100
    IDS 4
  )
endforeach()

set(_tests
  Fft::I4DivR1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 2500
    IDS 2
  )
endforeach()

set(_tests
  Fft::I2DivR1D1R1_ZeroTest:integrator
  Fft::I4DivR1D1R1_ZeroTest:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 200000
  )
endforeach()

set(_tests
  Fft::DivR1D1R1Test:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
  )
endforeach()

set(_tests
  Fft::DivR1Test:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS 4
  )
endforeach()

# Disable FFT tests temporarily
set_tests_properties(TransformWorlandTests_Fft_DivR1_integrator_id4 TransformWorlandTests_Fft_DivR1D1R1_integrator PROPERTIES DISABLED TRUE)
