message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
quicc_clear_test_list()
include(SetUlp)

quicc_set_ulp(_ulps
  ULP 40 40 40 70 600 800
  MPULP 40 40 40 60 600 800
  )

set(_tests
  P:<base_t>:projector
  DivY1:<base_t>:projector
  DivY2:<base_t>:projector
  P:<viewCpu_t>:projector
  DivY1:<viewCpu_t>:projector
  DivY2:<viewCpu_t>:projector
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 25 150 200 180 7000 70000
  MPULP 900 900 900 270 6000 70000
  )

set(_tests
  D1Y1:<base_t>:projector
  DivY1D1:<base_t>:projector
  DivY1D1Y1:<base_t>:projector
  DivY2D1Y1:<base_t>:projector
  D1Y1:<viewCpu_t>:projector
  D2Y1:<viewCpu_t>:projector
  DivY1D1:<viewCpu_t>:projector
  DivY1D1Y1:<viewCpu_t>:projector
  DivY1D2Y1:<viewCpu_t>:projector
  DivY2D1Y1:<viewCpu_t>:projector
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 10 160 450 1800 160000 5000000
  MPULP 10 150 650 2500 320000 5000000
  )

set(_tests
  D2Y1:<base_t>:projector
  DivY1D2Y1:<base_t>:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 40 200 625 1100 200000 9000000
  MPULP 600 600 600 25 25 25
  )

set(_tests
  SphRadLapl:<base_t>:projector
  SphRadLapl:<viewCpu_t>:projector
  )

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 1 40 70 180 12000 14000
  MPULP 150 165 525 25 25 25
  )

set(_tests
  D:<1,base_t>:projector
  D1:<viewCpu_t>:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 1 60 200 1800 85000 650000
  MPULP 150 165 525 1 1 1
  )

set(_tests
  D:<2,base_t>:projector
  D2:<viewCpu_t>:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 1 100 525 1200 12000000 25000000
  MPULP 150 165 525 1 1 1
  )

set(_tests
  D:<3,base_t>:projector
  D3:<viewCpu_t>:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 1 150 150 2500 8500000 300000000
  MPULP 150 165 525 1 1 1
  )

set(_tests
  D:<4,base_t>:projector
  D4:<viewCpu_t>:projector
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULPS ${_ulps}
    IDS 0 1 2 6 7 8
  )
endforeach()

# Merged tests
quicc_add_merged_test(
  ${TestExe}_projector_merged
  COMMAND ${TestExe}
  JN 4
)
