message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::PTest:projector
  Poly::DivR1Test:projector
  Poly::D1Test:projector
  Poly::DivR1D1R1Test:projector
  Poly::SphLaplTest:projector
  Poly::D1R1Test:projector
  Poly::CylLaplhTest:projector
  Poly::D1CylLaplhTest:projector
  Poly::DivR1_ZeroTest:projector
  Poly::D1_PTest:projector
  Poly::DivR1D1R1_ZeroTest:projector
  Poly::CylLaplh_DivR1D1R1Test:projector
  Poly::D1CylLaplh_D1DivR1D1R1Test:projector
  Poly::DivR1CylLaplh_ZeroTest:projector
  Fft::P:projector
  Fft::DivR1_Zero:projector
  Fft::D1:projector
  Fft::DivR1D1R1_Zero:projector
  Fft::SphLapl:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
  )
endforeach()
