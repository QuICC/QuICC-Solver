message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)

set(_tests
  Poly::D1:projector
  Poly::Ll:projector
  Poly::LlDivS1:projector
  Poly::LlDivS1Dp:projector
  Poly::LlD1:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 125
  )
endforeach()

set(_tests_templated
  Poly::P:<>:projector
  )

if(QUICC_USE_KOKKOS_CUDA)
list(APPEND _tests_templated
    Poly::PD1:<CudaIALegendreOperatorTypes>:projector
    Poly::PLl:<CudaIALegendreOperatorTypes>:projector
    Poly::PLlD1:<CudaIALegendreOperatorTypes>:projector
    Poly::PLlDivS1:<CudaIALegendreOperatorTypes>:projector
    Poly::PLlDivS1Dp:<CudaIALegendreOperatorTypes>:projector
  )
endif(QUICC_USE_KOKKOS_CUDA)

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 100
  )
endforeach()

set(_tests
  Poly::DivS1:projector
  Poly::DivS1Dp:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 261
  )
endforeach()

if(QUICC_USE_KOKKOS_CUDA)
set(_tests_templated
    Poly::PDivS1:<CudaIALegendreOperatorTypes>:projector
    Poly::PDivS1Dp:<CudaIALegendreOperatorTypes>:projector
  )

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 200
  )
endforeach()
endif(QUICC_USE_KOKKOS_CUDA)

# Special projector tests
set(_tests_templated
  Poly::P:<>:projector
  )

if(QUICC_USE_KOKKOS)
list(APPEND _tests_templated
  Poly::PP:<PIALegendreOperatorTypes>:projector
  )
endif(QUICC_USE_KOKKOS)

if(QUICC_USE_KOKKOS_CUDA)
list(APPEND _tests_templated
    Poly::PP:<CudaIALegendreOperatorTypes>:projector
  )
endif(QUICC_USE_KOKKOS_CUDA)

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  204 205  206   207
    ULPS 400 4500 19000 180000
    STEPS -1 -1 500 500
  )
endforeach()

# Perf test
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  108 108 108 108 108 108
    ULPS 20800 20800 20800 20800 20800 20800
    # node:rank, node*12:rank, node*36:rank
    SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
    STEPS 500 500 500 500 500 500
    PERFONLY
  )
endforeach()

# The following tests require ulps O(10^15) ...
# if(QUICC_MULTPRECISION)
  set(_disabled True)
# else()
  # set(_disabled False)
# endif()

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  208 209 210
    ULPS 150 150 500
    DISABLED ${_disabled}
  )
endforeach()

