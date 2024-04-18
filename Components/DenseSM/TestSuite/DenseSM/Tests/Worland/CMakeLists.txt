message(VERBOSE "Enabling Worland tests:")

set(TestExe DenseSMWorlandTests)

# Add target for all tests
add_executable(${TestExe} WorlandTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set(QUICC_WORK_DIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite")

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${QUICC_WORK_DIR}"
)

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT DenseSM
  FILENAME "Worland.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Create Geostrophic2Tor
set(_tests
  Geostrophic2Tor:dense
  )
quicc_set_ulp(_ulp
  ULP 50 70 11 18
  MPULP 25 25 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 100 101
    ULPS ${_ulp}
  )
endforeach()

# Create tests for GeostrophicAngular
set(_tests
  GeostrophicAngularMomentum:dense
  )
quicc_set_ulp(_ulp
  ULP 11 11 11 11 11
  MPULP 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4
    ULPS ${_ulp}
  )
endforeach()

# Create Tor2Geostrophic
set(_tests
  Tor2Geostrophic:dense
  )
quicc_set_ulp(_ulp
  ULP 11 35 11 40
  MPULP 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 100 101
    ULPS ${_ulp}
  )
endforeach()

# Create Tor2EGeostrophic
set(_tests
  Tor2EGeostrophic:dense
  )
quicc_set_ulp(_ulp
  ULP 11 11 11 11
  MPULP 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 100 101
    ULPS ${_ulp}
    DISABLED ON
  )
endforeach()

# Create Tor2GridS
set(_tests
  Tor2GridS:dense
  )
quicc_set_ulp(_ulp
  ULP 16 80 11 80
  MPULP 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 100 101
    ULPS ${_ulp}
  )
endforeach()

# Create Tor2Weights
set(_tests
  Tor2Weights:dense
  )
quicc_set_ulp(_ulp
  ULP 11 11
  MPULP 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
    DISABLED ON
  )
endforeach()

# Create tests for ProjFitEnergyR2 with Chebyshev type basis
set(_tests
  ProjFitEnergyR2:dense
  I2ProjFitEnergyR2:dense
  )
quicc_set_ulp(_ulp
  ULP 5000 2200 3800 1600 4500 1700 5900 4300 3200 2700 3000 1400
  MPULP 11 11 11 11 11 11 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4 5 6 7 8 9 10 11
    ULPS ${_ulp}
  )
endforeach()

# Create tests for ProjFitEnergyR2 with SphEnergy type basis
set(_tests
  ProjFitEnergyR2:dense
  I2ProjFitEnergyR2:dense
  )
quicc_set_ulp(_ulp
  ULP 5000 2200 3800 1600 4500 1700 5900 4300 3200 2700 3000 1400
  MPULP 11 11 11 11 11 11 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104 105 106 107 108 109 110 111
    ULPS ${_ulp}
  )
endforeach()

# Create tests for ProjFitEnergyPol with Chebyshev type basis
set(_tests
  ProjFitEnergyPol:dense
  I2ProjFitEnergyPol:dense
  )
quicc_set_ulp(_ulp
  ULP 9400 9800 6100 4500 3400 3600
  MPULP 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4 5
    ULPS ${_ulp}
  )
endforeach()

# Create tests for ProjFitEnergyPol with SphEnergy type basis
set(_tests
  ProjFitEnergyPol:dense
  I2ProjFitEnergyPol:dense
  )
quicc_set_ulp(_ulp
  ULP 9400 9800 6100 4500 3400 3600
  MPULP 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104 105
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
