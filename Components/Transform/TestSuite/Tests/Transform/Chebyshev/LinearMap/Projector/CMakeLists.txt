message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
include(SetUlp)

quicc_set_ulp(_ulps
  ULP 40 40 40 250 250 625
  MPULP 40 40 40 900 900 600
  )

set(_tests
  P:<base_t>:projector
  DivY1:<base_t>:projector
  DivY2:<base_t>:projector
  D1Y1:<base_t>:projector
  DivY1D1Y1:<base_t>:projector
  SphRadLapl:<base_t>:projector
  P:<viewCpu_t>:projector
  )

foreach(_test _ulp IN ZIP_LISTS _tests _ulps)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP ${_ulp}
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 60 375 425 116
  MPULP 150 165 525 150
  )

set(_tests
  D:<1,base_t>:projector
  D:<2,base_t>:projector
  D:<3,base_t>:projector
  D:<4,base_t>:projector
  )
foreach(_test _ulp IN ZIP_LISTS _tests _ulps)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Impl:Type
    ULP ${_ulp}
  )
endforeach()
