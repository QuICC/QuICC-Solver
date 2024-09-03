message(VERBOSE "Reductors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Value::Energy:<base_t>:reductor
  Value::EnergyR2:<base_t>:reductor
  Value::EnergySLaplR2:<base_t>:reductor
  Insulating::Energy:<base_t>:reductor
  Insulating::EnergyR2:<base_t>:reductor
  Insulating::EnergySLaplR2:<base_t>:reductor
  NoSlip::Energy:<base_t>:reductor
  NoSlip::EnergyR2:<base_t>:reductor
  NoSlip::EnergySLaplR2:<base_t>:reductor
)
quicc_set_ulp(_ulp_uni
  ULP 30 36 30 32
  MPULP 11 11 11 11
)
set(_ids_uni
  0 1 2 3
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

#set(_tests
#  Value::RadialPower:<base_t>:reductor
#  Value::RadialPowerDivR1:<base_t>:reductor
#  Value::RadialPowerDivR1D1R1:<base_t>:reductor
#  Insulating::RadialPower:<base_t>:reductor
#  Insulating::RadialPowerDivR1:<base_t>:reductor
#  Insulating::RadialPowerDivR1D1R1:<base_t>:reductor
#  NoSlip::RadialPower:<base_t>:reductor
#  NoSlip::RadialPowerDivR1:<base_t>:reductor
#  NoSlip::RadialPowerDivR1D1R1:<base_t>:reductor
#)
#quicc_set_ulp(_ulp_uni
#  ULP 60 120 120 120
#  MPULP 11 11 11 11
#)
#foreach(_test ${_tests})
#  quicc_add_test(${_test}
#    COMMAND ${TestExe}
#    KEYWORD Boundary::Operator:Template:type
#    IDS ${_ids_uni}
#    ULPS ${_ulp_uni}
#  )
#endforeach()

set(_tests
  Value::EnergyD1R1:<base_t>:reductor
  Insulating::EnergyD1R1:<base_t>:reductor
  NoSlip::EnergyD1R1:<base_t>:reductor
)
quicc_set_ulp(_ulp_uni
  ULP 30 201 230 240
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

set(_tests
  Value::Power:<base_t>:reductor
  Value::PowerR2:<base_t>:reductor
  Value::PowerD1R1:<base_t>:reductor
  Value::PowerSLaplR2:<base_t>:reductor
  Insulating::Power:<base_t>:reductor
  Insulating::PowerR2:<base_t>:reductor
  Insulating::PowerD1R1:<base_t>:reductor
  Insulating::PowerSLaplR2:<base_t>:reductor
  NoSlip::Power:<base_t>:reductor
  NoSlip::PowerR2:<base_t>:reductor
  NoSlip::PowerD1R1:<base_t>:reductor
  NoSlip::PowerSLaplR2:<base_t>:reductor
)
quicc_set_ulp(_ulp_uni
  ULP 30 500 1600 190
  MPULP 11 15 70 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests
  Value::EnergyR2:<base_t>:reductor
  Value::Energy:<base_t>:reductor
  Value::EnergySLaplR2:<base_t>:reductor
  Insulating::EnergyR2:<base_t>:reductor
  Insulating::Energy:<base_t>:reductor
  Insulating::EnergySLaplR2:<base_t>:reductor
  NoSlip::EnergyR2:<base_t>:reductor
  NoSlip::Energy:<base_t>:reductor
  NoSlip::EnergySLaplR2:<base_t>:reductor
)
set(_ids
  104 105 106 107
)
quicc_set_ulp(_ulp
  ULP 45 190 750 4750
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

#set(_tests
#  Value::RadialPower:<base_t>:reductor
#  Value::RadialPowerDivR1:<base_t>:reductor
#  Value::RadialPowerDivR1D1R1:<base_t>:reductor
#  Insulating::RadialPower:<base_t>:reductor
#  Insulating::RadialPowerDivR1:<base_t>:reductor
#  Insulating::RadialPowerDivR1D1R1:<base_t>:reductor
#  NoSlip::RadialPower:<base_t>:reductor
#  NoSlip::RadialPowerDivR1:<base_t>:reductor
#  NoSlip::RadialPowerDivR1D1R1:<base_t>:reductor
#)
#quicc_set_ulp(_ulp
#  ULP 280 900 3900 4850
#  MPULP 11 11 11 11
#)
#foreach(_test ${_tests})
#  quicc_add_test(${_test}
#    COMMAND ${TestExe}
#    KEYWORD Boundary::Operator:Template:type
#    IDS ${_ids}
#    ULPS ${_ulp}
#  )
#endforeach()

set(_tests
  Value::EnergyD1R1:<base_t>:reductor
  Insulating::EnergyD1R1:<base_t>:reductor
  NoSlip::EnergyD1R1:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 51 360 3305 10100
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Value::PowerD1R1:<base_t>:reductor
  Insulating::PowerD1R1:<base_t>:reductor
  NoSlip::PowerD1R1:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 160 350 525 9400
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Value::PowerR2:<base_t>:reductor
  Insulating::PowerR2:<base_t>:reductor
  NoSlip::PowerR2:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 125 280 880 2030
  MPULP 11 11 11 13
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Value::Power:<base_t>:reductor
  Insulating::Power:<base_t>:reductor
  NoSlip::Power:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 80 200 820 1970
  MPULP 11 11 11 21
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Value::PowerSLaplR2:<base_t>:reductor
  Insulating::PowerSLaplR2:<base_t>:reductor
  NoSlip::PowerSLaplR2:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 70 145 780 2140
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Boundary::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

# Split reductor tests

set(_split_tests
  Value::EnergyR2:<base_t>:reductor
  Insulating::EnergyR2:<base_t>:reductor
  NoSlip::EnergyR2:<base_t>:reductor
)

quicc_set_ulp(_ulp
  ULP 25 25
  MPULP 11 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Boundary::Operator:Template:type
        IDS  104 104
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 4:0 8:0
    )
endforeach()


quicc_set_ulp(_ulp
  ULP 3450 3200 3150 3050 3200 3450
  MPULP 11 11 11 11 11 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Boundary::Operator:Template:type
        IDS 108 108 108 108 108 108
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
        STEPS 500 500 500 500 500 1000
        PERFONLY
    )
endforeach()

set(_split_tests
  Value::PowerR2:<base_t>:reductor
  Insulating::PowerR2:<base_t>:reductor
  NoSlip::PowerR2:<base_t>:reductor
)

quicc_set_ulp(_ulp
  ULP 70 70 3750 3400 3350 3350 3350 3350
  MPULP 11 11 21 15 15 14 21 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Boundary::Operator:Template:type
        IDS  104 104 108 108 108 108 108 108
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
    )
endforeach()

#set(_split_tests
#  Value::RadialPower:<base_t>:reductor
#  Insulating::RadialPower:<base_t>:reductor
#  NoSlip::RadialPower:<base_t>:reductor
#)
#
#quicc_set_ulp(_ulp
#  ULP 280 280 35000 30500 30500 32500 30500 30500
#  MPULP 11 11 11 11 11 11 11 11
#)
#foreach(_test ${_split_tests})
#    quicc_add_test(${_test}
#        COMMAND ${TestExe}
#        KEYWORD Boundary::Operator:Template:type
#        IDS  104 104 108 108 108 108 108 108
#        ULPS ${_ulp}
#        # node:rank, node*12:rank, node*36:rank
#        SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
#    )
#endforeach()
