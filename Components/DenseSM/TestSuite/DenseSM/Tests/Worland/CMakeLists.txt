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
#quicc_fetch_test_reference(
#  ${TestExe}
#  COMPONENT DenseSM
#  FILENAME "Worland.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Create Geostrophic2Tor
set(_tests
  Geostrophic2Tor:dense
  PyGeostrophic2Tor:dense
  )
quicc_set_ulp(_ulp
  ULP 50 70
  MPULP 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1
    ULPS ${_ulp}
  )
endforeach()

# Create Geostrophic2Tor
set(_tests
  Geostrophic2Tor:dense
  PyGeostrophic2Tor:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500
  MPULP 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 2 3
    ULPS ${_ulp}
    DISABLED ON
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

# Create tests for ProjFitEnergy
set(_tests
  ProjFitEnergyPol:dense
  ProjFitEnergyR2:dense
  I2ProjFitEnergyR2:dense
  I2ProjFitEnergyPol:dense
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

list(POP_BACK CMAKE_MESSAGE_INDENT)
