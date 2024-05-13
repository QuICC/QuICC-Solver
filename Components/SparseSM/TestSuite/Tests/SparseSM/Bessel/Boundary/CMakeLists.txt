message(VERBOSE "Enabling Bessel boundary tests:")

set(TestExe SparseSMBesselBoundaryTests)

# Add target for all tests
add_executable(${TestExe} BesselBoundaryTests.cpp)

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
  Value:Value:boundary
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 0 1 2 3
    ULPS 100 100 100 100
  )
endforeach()

# Create basis tests targets
set(_tests
  Value:Insulating:boundary
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 0 1 2 3 4
    ULPS 100 100 100 100 100
  )
endforeach()

# Create basis tests targets
set(_tests
  D1:Value:boundary
  D1:Insulating:boundary
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 0 1 2 3 4
    ULPS 100 100 100 100 100
  )
endforeach()

# Create basis tests targets
set(_tests
  D2:Value:boundary
  D2:Insulating:boundary
  R1D1DivR1:Value:boundary
  R1D1DivR1:Insulating:boundary
  Value:Value:boundary
  Value:Insulating:boundary
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 0 1 2 3 4
    ULPS 120 120 120 280 400
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
