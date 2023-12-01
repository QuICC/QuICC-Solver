message(VERBOSE "Reductors:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

set(_tests
  Poly::Energy:<base_t>:reductor
  Poly::EnergyR2:<base_t>:reductor
  Poly::EnergySLaplR2:<base_t>:reductor
  Poly::Energy:<viewCpu_t>:reductor
  Poly::EnergyR2:<viewCpu_t>:reductor
  Poly::EnergySLaplR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::Energy:<viewGpu_t>:reductor
    Poly::EnergyR2:<viewGpu_t>:reductor
    Poly::EnergySLaplR2:<viewGpu_t>:reductor
  )
endif()
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
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

# triangular truncation
set(_tests
  Poly::Energy:<base_t>:reductor
  Poly::EnergyR2:<base_t>:reductor
  Poly::EnergySLaplR2:<base_t>:reductor
  # Poly::Energy:<viewCpu_t>:reductor
  # Poly::EnergyR2:<viewCpu_t>:reductor
  # Poly::EnergySLaplR2:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp_tri
  ULP 30 120 120 50
  MPULP 11 11 11 11
)
set(_ids_tri
  1000 1001 1002 1003
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_tri}
    ULPS ${_ulp_tri}
  )
endforeach()

set(_ids
  ${_ids_uni} ${_ids_tri}
)
set(_ulp
  ${_ulp_uni} ${_ulp_tri}
)

set(_tests
  Fft::EnergyR2:<base_t>:reductor
  Fft::Energy:<base_t>:reductor
  Fft::EnergySLaplR2:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::RadialPower:<base_t>:reductor
  Poly::RadialPowerDivR1:<base_t>:reductor
  Poly::RadialPowerDivR1D1R1:<base_t>:reductor
  Poly::RadialPower:<viewCpu_t>:reductor
  Poly::RadialPowerDivR1:<viewCpu_t>:reductor
  Poly::RadialPowerDivR1D1R1:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::RadialPower:<viewGpu_t>:reductor
    Poly::RadialPowerDivR1:<viewGpu_t>:reductor
    Poly::RadialPowerDivR1D1R1:<viewGpu_t>:reductor
  )
endif()
quicc_set_ulp(_ulp_uni
  ULP 60 120 120 120
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

set(_tests
  Poly::RadialPower:<base_t>:reductor
  Poly::RadialPowerDivR1:<base_t>:reductor
  Poly::RadialPowerDivR1D1R1:<base_t>:reductor
  # Poly::RadialPower:<viewCpu_t>:reductor
  # Poly::RadialPowerDivR1:<viewCpu_t>:reductor
  # Poly::RadialPowerDivR1D1R1:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp_tri
  ULP 60 220 120 130
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_tri}
    ULPS ${_ulp_tri}
  )
endforeach()

set(_tests
  Poly::EnergyD1R1:<base_t>:reductor
  Poly::EnergyD1R1:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp_uni
  ULP 30 201 230 240
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

# triangular truncation
set(_tests
  Poly::EnergyD1R1:<base_t>:reductor
  # Poly::EnergyD1R1:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp_tri
  ULP 30 200 230 240
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_tri}
    ULPS ${_ulp_tri}
  )
endforeach()

set(_ulp
  ${_ulp_uni} ${_ulp_tri}
)

set(_tests
  Fft::EnergyD1R1:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::Power:<base_t>:reductor
  Poly::PowerR2:<base_t>:reductor
  Poly::PowerD1R1:<base_t>:reductor
  Poly::PowerSLaplR2:<base_t>:reductor
  Poly::Power:<viewCpu_t>:reductor
  Poly::PowerR2:<viewCpu_t>:reductor
  Poly::PowerD1R1:<viewCpu_t>:reductor
  Poly::PowerSLaplR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::Power:<viewGpu_t>:reductor
    Poly::PowerR2:<viewGpu_t>:reductor
    Poly::PowerD1R1:<viewGpu_t>:reductor
    Poly::PowerSLaplR2:<viewGpu_t>:reductor
  )
endif()
quicc_set_ulp(_ulp_uni
  ULP 30 500 1600 190
  MPULP 11 15 70 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_uni}
    ULPS ${_ulp_uni}
  )
endforeach()

# triangular truncation
set(_tests
  Poly::PowerD1R1:<base_t>:reductor
  Poly::PowerR2:<base_t>:reductor
  Poly::Power:<base_t>:reductor
  Poly::PowerSLaplR2:<base_t>:reductor
  # Poly::PowerD1R1:<viewCpu_t>:reductor
  # Poly::PowerR2:<viewCpu_t>:reductor
  # Poly::Power:<viewCpu_t>:reductor
  # Poly::PowerSLaplR2:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp_tri
  ULP 30 500 1600 120
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids_tri}
    ULPS ${_ulp_tri}
  )
endforeach()


set(_ulp
  ${_ulp_uni} ${_ulp_tri}
)

set(_tests
  Fft::PowerD1R1:<base_t>:reductor
  Fft::PowerR2:<base_t>:reductor
  Fft::Power:<base_t>:reductor
  Fft::PowerSLaplR2:<base_t>:reductor
)

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests
  Poly::EnergyR2:<base_t>:reductor
  Poly::Energy:<base_t>:reductor
  Poly::EnergySLaplR2:<base_t>:reductor
  Poly::EnergyR2:<viewCpu_t>:reductor
  Poly::Energy:<viewCpu_t>:reductor
  Poly::EnergySLaplR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::Energy:<viewGpu_t>:reductor
    Poly::EnergyR2:<viewGpu_t>:reductor
    Poly::EnergySLaplR2:<viewGpu_t>:reductor
  )
endif()
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
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::EnergyR2:<base_t>:reductor
  Fft::Energy:<base_t>:reductor
  Fft::EnergySLaplR2:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::RadialPower:<base_t>:reductor
  Poly::RadialPowerDivR1:<base_t>:reductor
  Poly::RadialPowerDivR1D1R1:<base_t>:reductor
  Poly::RadialPower:<viewCpu_t>:reductor
  Poly::RadialPowerDivR1:<viewCpu_t>:reductor
  Poly::RadialPowerDivR1D1R1:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::RadialPower:<viewGpu_t>:reductor
    Poly::RadialPowerDivR1:<viewGpu_t>:reductor
    Poly::RadialPowerDivR1D1R1:<viewGpu_t>:reductor
  )
endif()
quicc_set_ulp(_ulp
  ULP 280 900 3900 4850
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::EnergyD1R1:<base_t>:reductor
  Poly::EnergyD1R1:<viewCpu_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 51 360 3305 10100
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::EnergyD1R1:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerD1R1:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 160 350 525 9400
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::PowerD1R1:<base_t>:reductor
)

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerR2:<base_t>:reductor
  Poly::PowerR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::PowerR2:<viewGpu_t>:reductor
  )
endif()
quicc_set_ulp(_ulp
  ULP 125 280 880 2030
  MPULP 11 11 11 13
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::PowerR2:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::Power:<base_t>:reductor
  Poly::Power:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _tests
    Poly::Power:<viewGpu_t>:reductor
  )
endif()
quicc_set_ulp(_ulp
  ULP 80 200 820 1970
  MPULP 11 11 11 21
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::Power:<base_t>:reductor
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerSLaplR2:<base_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 70 145 780 2140
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Fft::PowerSLaplR2:<base_t>:reductor
)

foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

# Kokkos reductor tests

if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
set(_tests
    Poly::EnergyR2:<kokkos_t>:reductor
    Poly::Energy:<kokkos_t>:reductor
    Poly::EnergySLaplR2:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 36 30 32
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::EnergyD1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 201 230 240
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerD1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 500 800 140
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerR2:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 300 1600 140
  MPULP 11 15 70 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
    Poly::RadialPower:<kokkos_t>:reductor
    Poly::RadialPowerDivR1:<kokkos_t>:reductor
    Poly::RadialPowerDivR1D1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 60 120 120 120
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
    Poly::Power:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 350 500 100
  MPULP 11 30 50 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerSLaplR2:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 30 140 250 65
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS 0 1 2 3
    ULPS ${_ulp}
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests
  Poly::EnergyR2:<kokkos_t>:reductor
  Poly::Energy:<kokkos_t>:reductor
  Poly::EnergySLaplR2:<kokkos_t>:reductor
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
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::RadialPower:<kokkos_t>:reductor
  Poly::RadialPowerDivR1:<kokkos_t>:reductor
  Poly::RadialPowerDivR1D1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 280 900 3900 4850
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::EnergyD1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 51 360 3305 10100
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerD1R1:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 160 350 525 9400
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerR2:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 125 280 880 2030
  MPULP 11 11 11 13
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::Power:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 80 200 820 1970
  MPULP 11 11 11 21
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

set(_tests
  Poly::PowerSLaplR2:<kokkos_t>:reductor
)
quicc_set_ulp(_ulp
  ULP 70 145 780 2140
  MPULP 11 11 11 11
)
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:Template:type
    IDS ${_ids}
    ULPS ${_ulp}
  )
endforeach()

endif()

# Split reductor tests

set(_split_tests
  Poly::EnergyR2:<base_t>:reductor
  Poly::EnergyR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _split_tests
    Poly::EnergyR2:<viewGpu_t>:reductor
  )
endif()
if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _split_tests
        Poly::EnergyR2:<kokkos_t>:reductor
        )
endif()

quicc_set_ulp(_ulp
  ULP 25 25
  MPULP 11 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
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
        KEYWORD Algorithm::Operator:Template:type
        IDS 108 108 108 108 108 108
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 8:0 96:0 288:0 256:0 3072:0 9216:0
        STEPS 500 500 500 500 500 1000
        PERFONLY
    )
endforeach()

set(_split_tests
  Poly::PowerR2:<base_t>:reductor
  Poly::PowerR2:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _split_tests
    Poly::PowerR2:<viewGpu_t>:reductor
  )
endif()

if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _split_tests
        Poly::PowerR2:<kokkos_t>:reductor
        )
endif()

quicc_set_ulp(_ulp
  ULP 70 70 3750 3400 3350 3350 3350 3350
  MPULP 11 11 21 15 15 14 21 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        IDS  104 104 108 108 108 108 108 108
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
    )
endforeach()

set(_split_tests
  Poly::RadialPower:<base_t>:reductor
  Poly::RadialPower:<viewCpu_t>:reductor
)
if(TARGET QuICC::Cuda)
  list(APPEND _split_tests
    Poly::RadialPower:<viewGpu_t>:reductor
  )
endif()
if(QUICC_USE_KOKKOS_CUDA OR QUICC_USE_KOKKOS_HIP)
    list(APPEND _split_tests
        Poly::RadialPower:<kokkos_t>:reductor
        )
endif()


quicc_set_ulp(_ulp
  ULP 280 280 35000 30500 30500 32500 30500 30500
  MPULP 11 11 11 11 11 11 11 11
)
foreach(_test ${_split_tests})
    quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Algorithm::Operator:Template:type
        IDS  104 104 108 108 108 108 108 108
        ULPS ${_ulp}
        # node:rank, node*12:rank, node*36:rank
        SPLITS 4:0 8:0 8:0 96:0 288:0 256:0 3072:0 9216:0
    )
endforeach()
