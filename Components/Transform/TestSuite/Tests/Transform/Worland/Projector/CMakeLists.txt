message(VERBOSE "Projectors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create  tests targets
include(ConfigureTests)

set(_tests
  Poly::P:<base_t>:projector
  Poly::P:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 8 94 55 60
  MPULP 1 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

# triangular truncation
set(_tests
  Poly::P:<base_t>:projector
  # Poly::P:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 8 140 55 60
  MPULP 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003
  )
endforeach()

# \todo DivR1 test missing

set(_tests
  Poly::D1:<base_t>:projector
  Poly::D1R1:<base_t>:projector
  Poly::DivR1_Zero:<base_t>:projector
  Poly::DivR1D1R1:<base_t>:projector
  Poly::DivR1D1R1_Zero:<base_t>:projector
  Poly::SphLapl:<base_t>:projector
  Poly::D1:<viewCpu_t>:projector
  Poly::D1R1:<viewCpu_t>:projector
  Poly::DivR1_Zero:<viewCpu_t>:projector
  Poly::DivR1D1R1:<viewCpu_t>:projector
  Poly::DivR1D1R1_Zero:<viewCpu_t>:projector
  Poly::SphLapl:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 9 60 49 51
  MPULP 1 1 2 2
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3
  )
endforeach()

# triangular truncation
set(_tests
  Poly::D1:<base_t>:projector
  Poly::D1R1:<base_t>:projector
  Poly::DivR1_Zero:<base_t>:projector
  Poly::DivR1D1R1:<base_t>:projector
  Poly::DivR1D1R1_Zero:<base_t>:projector
  Poly::SphLapl:<base_t>:projector
  # Poly::D1:<viewCpu_t>:projector
  # Poly::D1R1:<viewCpu_t>:projector
  # Poly::DivR1_Zero:<viewCpu_t>:projector
  # Poly::DivR1D1R1:<viewCpu_t>:projector
  # Poly::DivR1D1R1_Zero:<viewCpu_t>:projector
  # Poly::SphLapl:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 9 60 49 51
  MPULP 1 1 2 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 1000 1001 1002 1003
  )
endforeach()

set(_tests
  Poly::CylLaplh:<base_t>:projector
  Poly::CylLaplh_DivR1D1R1:<base_t>:projector
  Poly::D1_P:<base_t>:projector
  Poly::D1CylLaplh:<base_t>:projector
  Poly::D1CylLaplh_D1DivR1D1R1:<base_t>:projector
  Poly::DivR1:<base_t>:projector
  Poly::DivR1CylLaplh_Zero:<base_t>:projector
  Poly::CylLaplh:<viewCpu_t>:projector
  Poly::CylLaplh_DivR1D1R1:<viewCpu_t>:projector
  Poly::D1_P:<viewCpu_t>:projector
  Poly::D1CylLaplh:<viewCpu_t>:projector
  Poly::D1CylLaplh_D1DivR1D1R1:<viewCpu_t>:projector
  Poly::DivR1:<viewCpu_t>:projector
  Poly::DivR1CylLaplh_Zero:<viewCpu_t>:projector
)
quicc_set_ulp(_ulp
  ULP 1 1 1 1 1 1 1 1
  MPULP 1 1 1 1 1 1 1 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS ${_ulp}
    IDS 0 1 2 3 1000 1001 1002 1003
    DISABLED ON
  )
endforeach()

######################
# FFT based algorithm
######################

set(_tests
  Fft::P:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 4 17 5 5
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1_Zero:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 6 5 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::D1:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 4 4 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::DivR1D1R1_Zero:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 5 4 3 4
    IDS  0 1 2 3
  )
endforeach()

set(_tests
  Fft::SphLapl:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULPS 1 3 5 3
    IDS  0 1 2 3
  )
endforeach()

# Split projector tests

set(_split_tests_templated
  Poly::P:<base_t>:projector
  Poly::P:<viewCpu_t>:projector
)

foreach(_test ${_split_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  104 104 108 108 108 108 108 108
    ULPS 150 150 15000 15000 15000 15000 15000 15000
    # node:rank, node*12:rank, node*36:rank
    SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
  )
endforeach()
