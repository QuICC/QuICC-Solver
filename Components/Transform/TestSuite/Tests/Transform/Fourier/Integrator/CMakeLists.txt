message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)

# Mixed Fourier tests

set(_tests
  Mixed::P:<base_t>:integrator
  Mixed::D1:<base_t>:integrator
  Mixed::D1_P:<base_t>:integrator
  Mixed::D1_Neg:<base_t>:integrator
  Mixed::D2:<base_t>:integrator
  Mixed::P:<viewCpu_t>:integrator
  Mixed::D1:<viewCpu_t>:integrator
  Mixed::D1_P:<viewCpu_t>:integrator
  Mixed::D1_Neg:<viewCpu_t>:integrator
  Mixed::D2:<viewCpu_t>:integrator
)

if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Mixed::P:<viewGpu_t>:integrator
    Mixed::D1:<viewGpu_t>:integrator
    Mixed::D1_P:<viewGpu_t>:integrator
    Mixed::D1_Neg:<viewGpu_t>:integrator
    Mixed::D2:<viewGpu_t>:integrator
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Mixed::P:<viewGpuVkFFT_t>:integrator
    Mixed::D1:<viewGpuVkFFT_t>:integrator
    Mixed::D1_P:<viewGpuVkFFT_t>:integrator
    Mixed::D1_Neg:<viewGpuVkFFT_t>:integrator
    Mixed::D2:<viewGpuVkFFT_t>:integrator
  )
endif()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3 4 5 6
    ULP 70
  )
endforeach()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS   107 108  109
    ULPS  3600 7420 126000
  )
endforeach()

# Perf test
set(_perf_tests
  Mixed::P:<base_t>:integrator
  Mixed::D1:<base_t>:integrator
  Mixed::P:<viewCpu_t>:integrator
  Mixed::D1:<viewCpu_t>:integrator
)

if(QUICC_USE_CUFFT)
  list(APPEND _perf_tests
    Mixed::P:<viewGpu_t>:integrator
    Mixed::D1:<viewGpu_t>:integrator
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _perf_tests
    Mixed::P:<viewGpuVkFFT_t>:integrator
    Mixed::D1:<viewGpuVkFFT_t>:integrator
  )
endif()

foreach(_perf_id 108 109)
  if(${_perf_id} EQUAL 109)
    set(_perf_step_slow 50)
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
    list(APPEND _perf_ulps 7420)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_slow})
  endforeach()
  foreach(_s ${_perf_splits_fast})
    list(APPEND _perf_ids ${_perf_id})
    list(APPEND _perf_ulps 7420)
    list(APPEND _perf_splits ${_s})
    list(APPEND _perf_steps ${_perf_step_fast})
  endforeach()

  foreach(_test ${_perf_tests})
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Space::Operator:Template:type
      IDS ${_perf_ids}
      ULPS ${_perf_ulps}
      SPLITS ${_perf_splits}
      STEPS ${_perf_steps}
      PERFONLY
    )
  endforeach()
endforeach()

# Complex Fourier tests

set(_tests
  Complex::P:<base_t>:integrator
  Complex::P:<viewCpu_t>:integrator
  Complex::P_Clean:<base_t>:integrator
  Complex::P_Clean:<viewCpu_t>:integrator
  Complex::D1:<base_t>:integrator
  Complex::D1:<viewCpu_t>:integrator
  Complex::D1_P:<base_t>:integrator
  Complex::D1_P:<viewCpu_t>:integrator
  Complex::D1_Neg:<base_t>:integrator
  Complex::D1_Neg:<viewCpu_t>:integrator
  Complex::D2:<base_t>:integrator
  Complex::D2:<viewCpu_t>:integrator
  Complex::InvLapl2D:<base_t>:integrator
  Complex::InvLapl2D:<viewCpu_t>:integrator
  Complex::Df1InvLapl2D:<base_t>:integrator
  Complex::Df1InvLapl2D:<viewCpu_t>:integrator
  Complex::Mean:<base_t>:integrator
  Complex::Mean:<viewCpu_t>:integrator
)

if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Complex::P:<viewGpu_t>:integrator
    Complex::P_Clean:<viewGpu_t>:integrator
    Complex::D1:<viewGpu_t>:integrator
    Complex::D1_P:<viewGpu_t>:integrator
    Complex::D1_Neg:<viewGpu_t>:integrator
    Complex::D2:<viewGpu_t>:integrator
    Complex::InvLapl2D:<viewGpu_t>:integrator
    Complex::Df1InvLapl2D:<viewGpu_t>:integrator
    Complex::Mean:<viewGpu_t>:integrator
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Complex::P:<viewGpuVkFFT_t>:integrator
    Complex::P_Clean:<viewGpuVkFFT_t>:integrator
    Complex::D1:<viewGpuVkFFT_t>:integrator
    Complex::D1_P:<viewGpuVkFFT_t>:integrator
    Complex::D1_Neg:<viewGpuVkFFT_t>:integrator
    Complex::D2:<viewGpuVkFFT_t>:integrator
    Complex::InvLapl2D:<viewGpuVkFFT_t>:integrator
    Complex::Df1InvLapl2D:<viewGpuVkFFT_t>:integrator
    Complex::Mean:<viewGpuVkFFT_t>:integrator
  )
endif()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3 4 5 6
    ULP 70
  )
endforeach()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS   107  108   109
    ULPS  5168 16425 120000
  )
endforeach()


set(_tests
  Complex::Lapl2D:<base_t>:integrator
  Complex::Lapl2D:<viewCpu_t>:integrator
 )

if(QUICC_USE_CUFFT)
 list(APPEND _tests
   Complex::Lapl2D:<viewGpu_t>:integrator
 )
endif()
if(QUICC_USE_VKFFT)
 list(APPEND _tests
   Complex::Lapl2D:<viewGpuVkFFT_t>:integrator
 )
endif()


foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    ULP 90
  )
endforeach()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS   107  108   109
    ULPS  25520 81616 600000
  )
endforeach()
