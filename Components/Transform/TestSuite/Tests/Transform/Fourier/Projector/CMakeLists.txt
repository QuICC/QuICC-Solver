message(VERBOSE "Projectors:")

# Mixed Fourier tests

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Mixed::P:<base_t>:projector
  Mixed::P:<viewCpu_t>:projector
  Mixed::D1:<base_t>:projector
  Mixed::D1:<viewCpu_t>:projector
  Mixed::D2:<base_t>:projector
  Mixed::D2:<viewCpu_t>:projector
  Mixed::D3:<base_t>:projector
  Mixed::D3:<viewCpu_t>:projector
)

if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Mixed::P:<viewGpu_t>:projector
    Mixed::D1:<viewGpu_t>:projector
    Mixed::D2:<viewGpu_t>:projector
    Mixed::D3:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Mixed::P:<viewGpuVkFFT_t>:projector
    Mixed::D1:<viewGpuVkFFT_t>:projector
    Mixed::D2:<viewGpuVkFFT_t>:projector
    Mixed::D3:<viewGpuVkFFT_t>:projector
  )
endif()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3 4
    ULPS 70 70 150 13000 28000
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

set(_tests
  Mixed::D4:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:type
    IDS 0 1 2 3 4
    ULPS 70 70 250 250 1000
  )
endforeach()

# Perf test
set(_perf_tests
  Mixed::P:<base_t>:projector
  Mixed::P:<viewCpu_t>:projector
  Mixed::D1:<base_t>:projector
  Mixed::D1:<viewCpu_t>:projector
)

if(QUICC_USE_CUFFT)
  list(APPEND _perf_tests
    Mixed::P:<viewGpu_t>:projector
    Mixed::D1:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _perf_tests
    Mixed::P:<viewGpuVkFFT_t>:projector
    Mixed::D1:<viewGpuVkFFT_t>:projector
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
  Complex::P:<base_t>:projector
  Complex::P:<viewCpu_t>:projector
  Complex::D1:<base_t>:projector
  Complex::D1:<viewCpu_t>:projector
  Complex::D2:<base_t>:projector
  Complex::D2:<viewCpu_t>:projector
  Complex::D3:<base_t>:projector
  Complex::D3:<viewCpu_t>:projector
  Complex::D4:<base_t>:projector
  Complex::D4:<viewCpu_t>:projector
  Complex::Mean:<base_t>:projector
  Complex::Mean:<viewCpu_t>:projector
)

if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Complex::P:<viewGpu_t>:projector
    Complex::D1:<viewGpu_t>:projector
    Complex::D2:<viewGpu_t>:projector
    Complex::D3:<viewGpu_t>:projector
    Complex::D4:<viewGpu_t>:projector
    Complex::Mean:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Complex::P:<viewGpuVkFFT_t>:projector
    Complex::D1:<viewGpuVkFFT_t>:projector
    Complex::D2:<viewGpuVkFFT_t>:projector
    Complex::D3:<viewGpuVkFFT_t>:projector
    Complex::D4:<viewGpuVkFFT_t>:projector
    Complex::Mean:<viewGpuVkFFT_t>:projector
  )
endif()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3 4
    ULPS 190 190 190 190 1050
  )
endforeach()

foreach(_test ${_tests})
  string(REGEX MATCH "base_t" MATCHED ${_test})
  if(MATCHED)
    set(_disabled True)
  else()
    set(_disabled False)
  endif()
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS   107  108  109
    ULPS  4020 4020 17000
    DISABLED ${_disabled}
  )
endforeach()

set(_tests
  Complex::Lapl2D:<base_t>:projector
  Complex::Lapl2D:<viewCpu_t>:projector
  Complex::Df1Lapl2D:<viewCpu_t>:projector
  Complex::Ds1Lapl2D:<base_t>:projector
  Complex::Ds1Lapl2D:<viewCpu_t>:projector
)

if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Complex::Lapl2D:<viewGpu_t>:projector
    Complex::Df1Lapl2D:<viewGpu_t>:projector
    Complex::Ds1Lapl2D:<viewGpu_t>:projector
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Complex::Lapl2D:<viewGpuVkFFT_t>:projector
    Complex::Df1Lapl2D:<viewGpuVkFFT_t>:projector
    Complex::Ds1Lapl2D:<viewGpuVkFFT_t>:projector
  )
endif()

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3
    ULP 110
  )
endforeach()

set(_tests
  Complex::Df1Lapl2D:<base_t>:projector
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2 3
    ULPS 70 70 70 1100
  )
endforeach()

foreach(_test ${_tests})
  string(REGEX MATCH "base_t" MATCHED ${_test})
  if(MATCHED)
    set(_disabled True)
  else()
    set(_disabled False)
  endif()
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS   107   108   109
    ULPS  18600 18600 431017
    DISABLED ${_disabled}
  )
endforeach()
