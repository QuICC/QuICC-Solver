message(VERBOSE "BFLoop:")

# Create BFLoop tests targets
include(ConfigureTests)

set(_tests
  Mixed::P:<base_t>:bfloop
  Mixed::P:<viewCpu_t>:bfloop
  )
if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Mixed::P:<viewGpu_t>:bfloop
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Mixed::P:<viewGpuVkFFT_t>:bfloop
  )
endif()
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 0 1 2
    ULPS 11 11 11
  )
endforeach()

# Create BFLoop tests targets
set(_tests
  Mixed::P:<base_t>:bfloop
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 3 4 5 6
    ULPS 11 11 11 11
  )
endforeach()

# Create BFLoop tests targets
set(_tests
  Mixed::P:<viewCpu_t>:bfloop
  )
if(QUICC_USE_CUFFT)
  list(APPEND _tests
    Mixed::P:<viewGpu_t>:bfloop
  )
endif()
if(QUICC_USE_VKFFT)
  list(APPEND _tests
    Mixed::P:<viewGpuVkFFT_t>:bfloop
  )
endif()
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Space::Operator:Template:type
    IDS 3 4 5 6
    ULPS 11 11 11 11
    DISABLED ON
  )
endforeach()
