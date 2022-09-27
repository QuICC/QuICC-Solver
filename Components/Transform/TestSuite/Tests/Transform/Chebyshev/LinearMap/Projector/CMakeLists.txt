message(VERBOSE "Projectors:")

# Create projector tests targets
include(ConfigureTests)
include(SetUlp)

quicc_set_ulp(_ulps
  ULP 40 40 40 250 250 400
  MPULP 40 40 40 900 900 600
  )

set(_tests
  P:projector
  DivY1:projector
  DivY2:projector
  D1Y1:projector
  DivY1D1Y1:projector
  SphRadLapl:projector
  )

foreach(_test _ulp IN ZIP_LISTS _tests _ulps)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Type
    ULP ${_ulp}
  )
endforeach()

quicc_set_ulp(_ulps
  ULP 60 375 425 116
  MPULP 150 165 525 150
  )

set(_tests
  D:1:projector
  D:2:projector
  D:3:projector
  D:4:projector
  )
foreach(_test _ulp IN ZIP_LISTS _tests _ulps)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Template:Type
    ULP ${_ulp}
  )
endforeach()
