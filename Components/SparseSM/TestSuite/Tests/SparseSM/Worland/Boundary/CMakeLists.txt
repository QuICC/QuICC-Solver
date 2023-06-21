message(VERBOSE "Enabling Worland boundary tests:")

set(TestExe SparseSMWorlandBoundaryTests)

# Add target for all tests
add_executable(${TestExe} WorlandBoundaryTests.cpp)

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

# Create basis tests targets
set(_tests
  Value:boundary
  D1:boundary
  D2:boundary
  R1D1DivR1:boundary
  InsulatingSphere:boundary
  DivR1D1R1:boundary
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4 10 11 12 13 14
    ULP 40 40 40 40 40 40 40 40 40 40
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
