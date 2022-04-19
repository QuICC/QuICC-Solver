message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)
#set(_tests
#  Poly::R1Test:integrator
#  Poly::PTest:integrator
#  Poly::DivR1Test:integrator
#  Poly::DivR1D1R1Test:integrator
#  Poly::I2_ZeroTest:integrator
#  Poly::I2DivR1_ZeroTest:integrator
#  Poly::I4DivR1_ZeroTest:integrator
#  Poly::I2DivR1D1R1_ZeroTest:integrator
#  Poly::I4DivR1D1R1_ZeroTest:integrator
#  Poly::I4DivR1D1R1_I2Test:integrator
#  Poly::I6CylLaplh_I4D1R1Test:integrator
#  Poly::I6DivR1D1R1_I4Test:integrator
#  Poly::I6DivR1_ZeroTest:integrator
#  )
set(_tests
  Poly::R1Test:integrator
  Poly::PTest:integrator
  Poly::DivR1Test:integrator
  Fft::R1Test:integrator
  Fft::PTest:integrator
  Fft::DivR1Test
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
  )
endforeach()
