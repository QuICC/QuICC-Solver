message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)

set(_tests_templated
    Poly::P:<base_t>:integrator
    Poly::D1:<base_t>:integrator
    Poly::DivLl:<base_t>:integrator
    Poly::DivLlD1:<base_t>:integrator
    Poly::DivLlDivS1:<base_t>:integrator
    Poly::DivLlDivS1Dp:<base_t>:integrator
    Poly::DivS1:<base_t>:integrator
    Poly::DivS1Dp:<base_t>:integrator
    Poly::Ll:<base_t>:integrator
    Poly::LlD1:<base_t>:integrator
    Poly::LlDivS1:<base_t>:integrator
    Poly::LlDivS1Dp:<base_t>:integrator
    Poly::P:<viewCpu_t>:integrator
    Poly::D1:<viewCpu_t>:integrator
    Poly::DivLl:<viewCpu_t>:integrator
    Poly::DivLlD1:<viewCpu_t>:integrator
    Poly::DivLlDivS1:<viewCpu_t>:integrator
    Poly::DivLlDivS1Dp:<viewCpu_t>:integrator
    Poly::DivS1:<viewCpu_t>:integrator
    Poly::DivS1Dp:<viewCpu_t>:integrator
    Poly::Ll:<viewCpu_t>:integrator
    Poly::LlD1:<viewCpu_t>:integrator
    Poly::LlDivS1:<viewCpu_t>:integrator
    Poly::LlDivS1Dp:<viewCpu_t>:integrator
)

if(TARGET QuICC::Cuda)
    list(APPEND _tests_templated
        Poly::P:<viewGpu_t>:integrator
        Poly::D1:<viewGpu_t>:integrator
        Poly::DivLl:<viewGpu_t>:integrator
        Poly::DivLlD1:<viewGpu_t>:integrator
        Poly::DivLlDivS1:<viewGpu_t>:integrator
        Poly::DivLlDivS1Dp:<viewGpu_t>:integrator
        Poly::DivS1:<viewGpu_t>:integrator
        Poly::DivS1Dp:<viewGpu_t>:integrator
        Poly::Ll:<viewGpu_t>:integrator
        Poly::LlD1:<viewGpu_t>:integrator
        Poly::LlDivS1:<viewGpu_t>:integrator
        Poly::LlDivS1Dp:<viewGpu_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _tests_templated
        Poly::DivLl:<kokkos_t>:integrator
        Poly::DivS1Dp:<kokkos_t>:integrator
        Poly::DivLlD1:<kokkos_t>:integrator
        Poly::DivS1:<kokkos_t>:integrator
        Poly::DivLlDivS1:<kokkos_t>:integrator
        Poly::DivLlDivS1Dp:<kokkos_t>:integrator
        Poly::LlDivS1:<kokkos_t>:integrator
        Poly::LlDivS1Dp:<kokkos_t>:integrator
        Poly::D1:<kokkos_t>:integrator
        Poly::Ll:<kokkos_t>:integrator
        Poly::LlD1:<kokkos_t>:integrator
    )
endif()

foreach(_test ${_tests_templated})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        ULP 150
)
endforeach()


set(_tests_templated
    Poly::Ll2:<base_t>:integrator
)

if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _tests_templated
        Poly::Ll2:<kokkos_t>:integrator
    )
endif()
foreach(_test ${_tests_templated})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        ULP 1600
)
endforeach()

# Special integrator tests

set(_tests_templated
    Poly::P:<base_t>:integrator
    Poly::P:<viewCpu_t>:integrator
)

if(TARGET QuICC::Cuda)
    list(APPEND _tests_templated
        Poly::P:<viewGpu_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS AND NOT QUICC_USE_KOKKOS_CUDA)
    list(APPEND _tests_templated
        Poly::P:<kokkos_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _tests_templated
        Poly::P:<kokkos_t>:integrator
    )
endif()

# Special integrator tests
foreach(_test ${_tests_templated})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        IDS  204 205 206 207  208
        ULPS 30  120 370 1612 7925
    )
endforeach()

# Split integrator tests

set(_split_tests_templated
    Poly::P:<base_t>:integrator
    Poly::P:<viewCpu_t>:integrator
)

if(TARGET QuICC::Cuda)
    list(APPEND _split_tests_templated
        Poly::P:<viewGpu_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS AND NOT QUICC_USE_KOKKOS_CUDA)
    list(APPEND _split_tests_templated
        Poly::P:<kokkos_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS_CUDA)
    list(APPEND _split_tests_templated
        Poly::P:<kokkos_t>:integrator
    )
endif()

foreach(_test ${_split_tests_templated})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        IDS  204 204 204 204 204 204
        ULPS 50 50 50 50 50 50
        SPLITS 4:0 4:1 4:3 8:0 8:3 8:7
    )
endforeach()

# Split integrator performance tests

foreach(_test ${_split_tests_templated})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        IDS  108 108 108 108 108 108
        ULPS 20800 20800 20800 20800 20800 20800
        # node:rank, node*12:rank, node*36:rank
        SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
        STEPS 500 500 500 500 500 1000
        PERFONLY
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
