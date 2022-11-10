message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
set(_tests
  Poly::P:projector
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

# Special projector tests
set(_tests
  Poly::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS  204 205  206   207
    ULPS 400 4500 19000 180000
    STEPS -1 -1 500 500
  )
endforeach()

# Perf test
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
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

# Special projector tests
set(_tests
  Poly::P:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    IDS  208 209 210
    ULPS 150 150 500
    DISABLED ${_disabled}
  )
endforeach()
