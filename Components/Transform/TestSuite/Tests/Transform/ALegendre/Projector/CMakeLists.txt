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
    ULP 200
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
