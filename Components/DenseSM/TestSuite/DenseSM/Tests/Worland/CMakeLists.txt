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
#
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Create Geostrophic2Tor
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR1D1F:dense
  DivR1D1CF:dense
  DivR1D1FC:dense
  DivR2D1R1F:dense
  DivR2CFD1R1:dense
  DivR2D1R1FC:dense
  DivR2FD1R1:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 550 750 500 1500
  MPULP 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 10 11
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
