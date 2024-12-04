message(VERBOSE "Enabling Chebyshev tests:")

set(TestExe DenseSMChebyshevTests)

# Add target for all tests
add_executable(${TestExe} ChebyshevTests.cpp)

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
#  FILENAME "Chebyshev.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)
#
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Create Geostrophic2Tor
set(_tests
  R2DivR1F:dense
  R4DivR1F:dense
  R4DivR1D1F:dense
  R2DivR2D1R1F:dense
  R4DivR2D1R1F:dense
  R2DivR2FD1R1:dense
  R4DivR2FD1R1:dense
  R4DivR1D1DivR1FD1R1:dense
  R4DivR1D1DivR1D1R1F:dense
  R4DivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500
  MPULP 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
