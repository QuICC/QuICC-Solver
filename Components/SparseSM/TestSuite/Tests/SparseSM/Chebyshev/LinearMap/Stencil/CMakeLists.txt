message(VERBOSE "Enabling Chebyshev LinearMap boundary stencil tests:")

set(TestExe SparseSMChebyshevLinearMapStencilTests)

# Add target for all tests
add_executable(${TestExe} ChebyshevLinearMapStencilTests.cpp)

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
  Value:stencil
  D1:stencil
  R1D1DivR1:stencil
  InsulatingShell:stencil
  ValueD1:stencil
  ValueD2:stencil
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 10 11
    ULPS 11 11 11 11
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
