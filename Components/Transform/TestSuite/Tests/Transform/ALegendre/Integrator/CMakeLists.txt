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
set(_perf_tests
    Poly::P:<base_t>:integrator
    Poly::P:<viewCpu_t>:integrator
)

if(TARGET QuICC::Cuda)
    list(APPEND _perf_tests
        Poly::P:<viewGpu_t>:integrator
    )
endif()

if(QUICC_USE_KOKKOS)
    list(APPEND _perf_tests
        Poly::P:<kokkos_t>:integrator
    )
endif()

foreach(_perf_id 108 109)
  if(${_perf_id} EQUAL 109)
    set(_perf_step_slow 10)
    set(_perf_step_fast 50)
  else()
    set(_perf_step_slow 500)
    set(_perf_step_fast 1000)
  endif()

  set(_perf_splits_slow 8:0 16:0 32:0 64:0 96:0 128:0 192:0 256:0 288:0 384:0)
  set(_perf_splits_fast 576:0 768:0 1152:0 1536:0 2304:0 3072:0 4608:0 9216:0)

  set(_perf_ids )
  set(_perf_ulps )
  set(_perf_steps )
  set(_perf_splits )
  foreach(_s ${_perf_splits_slow})
    list(APPEND _perf_ids ${_perf_id})
    list(APPEND _perf_ulps 20800)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_slow})
  endforeach()
  foreach(_s ${_perf_splits_fast})
    list(APPEND _perf_ids ${_perf_id})
    list(APPEND _perf_ulps 20800)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_fast})
  endforeach()

  foreach(_test ${_perf_tests})
      quicc_add_test(${_test}
          COMMAND ${TestExe}
          KEYWORD Algorithm::Operator:Template:type
          IDS ${_perf_ids}
          ULPS ${_perf_ulps}
          SPLITS ${_perf_splits}
          STEPS ${_perf_steps}
          PERFONLY
      )
  endforeach()
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
