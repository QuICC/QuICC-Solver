message(VERBOSE "Enabling Chebyshev LinearMap tests:")

set(TestExe SparseSMChebyshevLinearMapTests)

# Add target for all tests
add_executable(${TestExe} ChebyshevLinearMapTests.cpp)

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
  COMPONENT SparseSM
  FILENAME "ChebyshevLinearMap.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Create basis tests targets
include(ConfigureTests)
set(_tests
  Y1:sparse
  Y2:sparse
  I1:sparse
  I2:sparse
  I2Y1:sparse
  I2Y2:sparse
  I2D1:sparse
  I2Y1D1Y1:sparse
  I3:sparse
  I4:sparse
  I4D1:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 10 11
    ULP 11
  )
endforeach()

# Create basis tests targets
set(_tests
  I4Y3:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 10 11
    ULPS 25 25 11 11
  )
endforeach()

# Create basis tests targets
set(_tests
  I4Y3D1Y1:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 10 11
    ULPS 18 58 11 36
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
