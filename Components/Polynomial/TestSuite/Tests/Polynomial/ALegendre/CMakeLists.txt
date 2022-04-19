message(VERBOSE "Enabling ALegendre tests:")

set(TestExe PolynomialALegendreTests)

# Add target for all tests
add_executable(${TestExe} ALegendreTests.cpp)

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

# Create ALegendre tests targets
include(ConfigureTests)
set(ALegendreTests PlmTest dPlmTest sin_1PlmTest)
foreach(_test ${ALegendreTests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator
    TYPES matrix weighted
  )
endforeach()

# Base dir for reference data
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Polynomial
  FILENAME "ALegendre.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)
