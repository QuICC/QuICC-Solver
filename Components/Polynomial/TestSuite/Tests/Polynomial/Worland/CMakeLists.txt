message(VERBOSE "Enabling Worland tests:")

set(TestExe PolynomialWorlandTests)

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

# Create Worland tests targets
include(ConfigureTests)
set(WorlandTests
  Wnl
  dWnl
  slaplWnl
  drWnl
  dr_1drWnl
  r_1claplhWnl
  claplhWnl
  dclaplhWnl
)
include(SetUlp)
quicc_set_ulp(_ulp
  ULP 16000
  MPULP 11
)
foreach(_test ${WorlandTests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator
    TYPES matrix weighted
    ULP ${_ulp}
  )
endforeach()

set(WorlandTestsTemplated
  r_1Wnl:recurrence_t
  r_1Wnl:explicit_t
)
quicc_set_ulp(_ulp
  ULP 10000
  MPULP 11
)
foreach(_test ${WorlandTestsTemplated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Template
    TYPES matrix weighted
    ULP ${_ulp}
  )
endforeach()

set(WorlandTestsTemplated
  r_1drWnl:recurrence_t
)
quicc_set_ulp(_ulp
  ULP 10000
  MPULP 11
)
foreach(_test ${WorlandTestsTemplated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Template
    TYPES matrix
    ULP ${_ulp}
  )
endforeach()

quicc_set_ulp(_ulp
  ULP 20000
  MPULP 11
)
set(WorlandTestsTemplated
  r_1Wnl:implicit_t
  r_1drWnl:explicit_t
  r_1drWnl:implicit_t
)
foreach(_test ${WorlandTestsTemplated})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:Template
    TYPES weighted
    ULP ${_ulp}
  )
endforeach()

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Polynomial
  FILENAME "Worland.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)
