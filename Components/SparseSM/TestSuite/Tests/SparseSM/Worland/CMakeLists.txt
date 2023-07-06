message(VERBOSE "Enabling Worland tests:")

set(TestExe SparseSMWorlandTests)

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
  COMPONENT SparseSM
  FILENAME "Worland.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Create Worland basis tests targets
include(ConfigureTests)
set(_tests
  R2:sparse
  I2:sparse
  I2Lapl:sparse
  I2Qp:sparse
  I4:sparse
  I4Lapl:sparse
  I4Lapl2:sparse
  I4Qp:sparse
  I6:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4
    ULP 70
  )
endforeach()
set(_tests
  I2Qm:sparse
  I4Qm:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3
    ULP 70
  )
endforeach()

# Create Worland basis tests targets
set(_tests
  R2:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 10 11 12 13 14
    ULP 70
  )
endforeach()

# Create Worland basis tests targets
set(_tests
  I2:sparse
  I2Lapl:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 10 11 12 13 14
    ULP 70
  )
endforeach()

# Create Worland basis tests targets
set(_tests
  I2Qm:sparse
  I4Qm:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 10 11 12 13
    ULP 75
  )
endforeach()

# Create Worland basis tests targets
set(_tests
  I4:sparse
  I4Lapl:sparse
  I4Lapl2:sparse
  I6:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 10 11 12 13 14
    ULP 95
  )
endforeach()

# Create tests for truncated operators
set(_tests
  I2:sparse
  I2Qp:sparse
  I4:sparse
  I4Lapl:sparse
  I4Qp:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104
    ULP 70
  )
endforeach()
set(_tests
  I4Qm:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103
    ULP 70
  )
endforeach()

# Create tests for truncated operators
set(_tests
  I2:sparse
  I2Qp:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 110 111 112 113 114
    ULP 80
  )
endforeach()

# Create tests for truncated operators
set(_tests
  I4:sparse
  I4Lapl:sparse
  I4Qp:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 110 111 112 113 114
    ULP 100
  )
endforeach()

set(_tests
  I4Qm:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 110 111 112 113
    ULP 75
  )
endforeach()

add_subdirectory(Boundary)

list(POP_BACK CMAKE_MESSAGE_INDENT)
