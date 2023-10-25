message(VERBOSE "Projectors:")

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

# Perf test
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS  108 108 108 108 108 108
    ULPS 7420 7420 7420 7420 7420 7420
    # node:rank, node*12:rank, node*36:rank
    SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
    STEPS   100   500   500   500    500   1000
    PERFONLY
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
