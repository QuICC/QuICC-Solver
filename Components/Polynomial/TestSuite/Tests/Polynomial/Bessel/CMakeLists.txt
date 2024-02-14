message(VERBOSE "Enabling Bessel tests:")

set(TestExe PolynomialBesselTests)

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

# Create Bessel tests targets
include(ConfigureTests)
set(BesselTests
  TorSphJnl
  PolSphJnl
  rTorSphJnl
  rPolSphJnl
  r_1TorSphJnl
  r_1PolSphJnl
  dTorSphJnl
  dPolSphJnl
  drTorSphJnl
  drPolSphJnl
  slaplTorSphJnl
  slaplPolSphJnl
)
include(SetUlp)
quicc_set_ulp(_ulp
  ULP 16000
  MPULP 11
)
foreach(_test ${BesselTests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator
    TYPES matrix weighted
    ULP ${_ulp}
  )
endforeach()

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
#quicc_fetch_test_reference(
#  ${TestExe}
#  COMPONENT Polynomial
#  FILENAME "Bessel.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)
