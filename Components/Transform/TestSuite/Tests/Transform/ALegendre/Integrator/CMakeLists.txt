message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)

set(_tests
  Poly::DivLlDivS1DpTest:integrator
  Poly::DivLlDivS1Test:integrator
  Poly::DivLl:integrator
  Poly::DivS1Dp:integrator
  Poly::DivS1:integrator
  Poly::LlDivS1Dp:integrator
  Poly::LlDivS1:integrator
  Poly::Ll:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 130
  )
endforeach()

set(_tests_templated
  Poly::P:<>:integrator
  )

if(QUICC_USE_KOKKOS_CUDA)
list(APPEND _tests_templated
    Poly::PDivLl:<CudaIALegendreOperatorTypes>:integrator
    Poly::PDivS1Dp:<CudaIALegendreOperatorTypes>:integrator
    Poly::PDivLlD1:<CudaIALegendreOperatorTypes>:integrator
    Poly::PDivS1:<CudaIALegendreOperatorTypes>:integrator
    Poly::PDivLlDivS1:<CudaIALegendreOperatorTypes>:integrator
    Poly::PDivLlDivS1Dp:<CudaIALegendreOperatorTypes>:integrator
    Poly::PLlDivS1:<CudaIALegendreOperatorTypes>:integrator
    Poly::PLlDivS1Dp:<CudaIALegendreOperatorTypes>:integrator
   )
endif(QUICC_USE_KOKKOS_CUDA)
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 150
  )
endforeach()

set(_tests
  Poly::D1Test:integrator
  Poly::LlD1:integrator
  Poly::Ll:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 150
  )
endforeach()

if(QUICC_USE_KOKKOS_CUDA)
set(_tests_templated
    Poly::PD1:<CudaIALegendreOperatorTypes>:integrator
    Poly::PLl:<CudaIALegendreOperatorTypes>:integrator
    Poly::PLlD1:<CudaIALegendreOperatorTypes>:integrator
  )

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 150
  )
endforeach()
endif(QUICC_USE_KOKKOS_CUDA)

set(_tests
  Poly::Ll2:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULP 1400
  )
endforeach()

if(QUICC_USE_KOKKOS_CUDA)
set(_tests_templated
    Poly::PLl2:<CudaIALegendreOperatorTypes>:integrator
   )
foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    ULP 1600
  )
endforeach()
endif(QUICC_USE_KOKKOS_CUDA)

# Special integrator tests
set(_tests_templated
  Poly::P:<>:integrator
  )
if(QUICC_USE_KOKKOS)
list(APPEND _tests_templated
  Poly::PP:<PIALegendreOperatorTypes>:integrator
  )
endif(QUICC_USE_KOKKOS)

if(QUICC_USE_KOKKOS_CUDA)
list(APPEND _tests_templated
    Poly::PP:<CudaIALegendreOperatorTypes>:integrator
  )
endif(QUICC_USE_KOKKOS_CUDA)

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  204 205 206 207  208
    ULPS 30  120 370 1612 7925
  )
endforeach()

# The following datasets need to be generated
# if(QUICC_MULTPRECISION)
  set(_disabled True)
# else()
#   set(_disabled False)
# endif()

foreach(_test ${_tests_templated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS  209    210
    ULPS 15000  80000
    DISABLED ${_disabled}
  )
endforeach()
