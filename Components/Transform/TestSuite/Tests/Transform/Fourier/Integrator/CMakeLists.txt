message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)
set(_tests
  Complex::P:integrator
  Complex::D1:integrator
  Complex::D1_Neg:integrator
  Complex::D1_P:integrator
  Complex::D2:integrator
  Complex::Lapl2D:integrator
  Complex::InvLapl2D:integrator
  Complex::Df1InvLapl2D:integrator
  Complex::Mean:integrator
  Complex::P_Clean:integrator
  Mixed::P:integrator
  Mixed::D1:integrator
  Mixed::D1_Neg:integrator
  Mixed::D1_P:integrator
  Mixed::D2:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
  )
endforeach()
