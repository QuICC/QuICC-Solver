message(VERBOSE "Enabling Jacobi tests:")

set(TestExe PolynomialJacobiTests)

# Add target for all tests
add_executable(${TestExe} JacobiTests.cpp)

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

# Create Jacobi tests targets
include(ConfigureTests)
set(JacobiTests PnabTest dPnabTest d2PnabTest d3PnabTest)
foreach(_test ${JacobiTests})
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
  FILENAME "Jacobi.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

add_executable(JacobiBaseTest JacobiBaseTests.cpp)
target_link_libraries(JacobiBaseTest ${QUICC_CURRENT_COMPONENT_LIB} Catch2::Catch2 )

add_test(
    NAME JacobiBaseTest
    COMMAND JacobiBaseTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )

if(NOT QUICC_MULTPRECISION)
  add_executable(JacobiAsymptoticTest JacobiAsymptoticTests.cpp)
  target_link_libraries(JacobiAsymptoticTest ${QUICC_CURRENT_COMPONENT_LIB} Catch2::Catch2 )

  add_test(
      NAME JacobiAsymptoticTest
      COMMAND JacobiAsymptoticTest
      WORKING_DIRECTORY
      "${QUICC_WORK_DIR}"
    )
endif()
