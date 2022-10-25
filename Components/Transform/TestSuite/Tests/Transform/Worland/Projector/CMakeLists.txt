message(VERBOSE "Projectors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create  tests targets
include(ConfigureTests)

set(_tests
  Poly::P:projector
  )
quicc_set_ulp(_ulp
  ULP 8 94 55 60
  MPULP 1 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::DivR1_Zero:projector
  )
quicc_set_ulp(_ulp
  ULP 8 48 49 51
  MPULP 1 1 2 2
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::D1:projector
  )
quicc_set_ulp(_ulp
  ULP 8 60 40 39
  MPULP 1 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::DivR1D1R1_Zero:projector
  )
quicc_set_ulp(_ulp
  ULP 9 60 40 42
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::SphLapl:projector
  )
quicc_set_ulp(_ulp
  ULP 1 35 32 31
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::DivR1D1R1:projector
  )
quicc_set_ulp(_ulp
  ULP 9 60 40 39
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::D1R1:projector
  )
quicc_set_ulp(_ulp
  ULP 9 54 42 41
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

set(_tests
  Poly::CylLaplh:projector
  Poly::CylLaplh_DivR1D1R1:projector
  Poly::D1CylLaplh:projector
  Poly::D1_P:projector
  Poly::D1CylLaplh_D1DivR1D1R1:projector
  Poly::DivR1CylLaplh_Zero:projector
  )
quicc_set_ulp(_ulp
  ULP 1 1 1 1
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3
    DISABLED ON
  )
endforeach()

######################
# FFT based algorithm
######################

set(_tests
  Fft::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 4 17 5 5
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1_Zero:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 5 6 5 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::D1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 5 4 4 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1D1R1_Zero:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 5 4 3 3
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::SphLapl:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 1 3 5 3
    IDS  0 1 2 3
  )
endforeach()
