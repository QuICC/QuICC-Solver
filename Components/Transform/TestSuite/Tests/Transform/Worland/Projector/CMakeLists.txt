message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::PTest:projector
  Poly::CylLaplhTest:projector
  Poly::CylLaplh_DivR1D1R1Test:projector
  Poly::DivR1CylLaplh_ZeroTest:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 5001
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::DivR1Test:projector
  Poly::DivR1_ZeroTest:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 3001
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::D1Test:projector
  Poly::DivR1D1R1Test:projector
  Poly::SphLaplTest:projector
  Poly::D1R1Test:projector
  Poly::D1CylLaplhTest:projector
  Poly::D1_PTest:projector
  Poly::DivR1D1R1_ZeroTest:projector
  Poly::D1CylLaplh_D1DivR1D1R1Test:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1000
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Fft::DivR1_Zero:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 750
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Fft::D1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1800
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Fft::SphLapl:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 16000
  )
endforeach()

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Fft::P:projector
  Fft::DivR1D1R1_Zero:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 3500
  )
endforeach()
