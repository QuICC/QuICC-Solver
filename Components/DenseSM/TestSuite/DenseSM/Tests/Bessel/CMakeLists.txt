message(VERBOSE "Enabling Bessel tests:")

set(TestExe DenseSMBesselTests)

# Add target for all tests
add_executable(${TestExe} BesselTests.cpp)

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
  FILENAME "Bessel.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Create tests for ProjFitEnergyR2
set(_tests
  ProjFitEnergyR2:dense
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

# Create tests for ProjFitEnergyPol
set(_tests
  ProjFitEnergyPol:dense
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

# Create tests for Coriolis
set(_tests
  CoriolisQm:dense
  )
quicc_set_ulp(_ulp
  ULP 25 25 25 25 25 25 25 25 25 25 35
  MPULP 11 11 11 11 11 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 1000 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010
    ULPS ${_ulp}
  )
endforeach()

# Create tests for Coriolis
set(_tests
  CoriolisQp:dense
  )
quicc_set_ulp(_ulp
  ULP 1600 25 600 25 250 25 450 25 320 25 110
  MPULP 11 11 11 11 11 11 11 11 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 1000 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010
    ULPS ${_ulp}
  )
endforeach()

# Create Geostrophic2Tor
set(_tests
  Geostrophic2Tor:dense
  )
quicc_set_ulp(_ulp
  ULP 20 60
  MPULP 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
  )
endforeach()

# Create GeostrophicAngularMomentum
set(_tests
  GeostrophicAngularMomentum:dense
  )
quicc_set_ulp(_ulp
  ULP 1500 11 30
  MPULP 11 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2
    ULPS ${_ulp}
  )
endforeach()

# Create GeostrophicEnergy
set(_tests
  GeostrophicEnergy:dense
  )
quicc_set_ulp(_ulp
  ULP 15 50
  MPULP 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
  )
endforeach()

# Create Tor2GridS
set(_tests
  Tor2GridS:dense
  )
quicc_set_ulp(_ulp
  ULP 25 400
  MPULP 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
  )
endforeach()

# Create Tor2Geostrophic
set(_tests
  Tor2Geostrophic:dense
  )
quicc_set_ulp(_ulp
  ULP 20 37
  MPULP 11 11
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
