message(VERBOSE "Integrators:")

# Create integrator tests targets
include(ConfigureTests)

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
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS  108 108 108 108 108 108
    ULPS 7420 7420 7420 7420 7420 7420
    # node:rank, node*12:rank, node*36:rank
    SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
    STEPS   25   50   100   100    500   1000
    PERFONLY
  )
endforeach()


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
