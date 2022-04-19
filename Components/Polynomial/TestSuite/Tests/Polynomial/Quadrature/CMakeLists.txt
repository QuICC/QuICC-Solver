message(VERBOSE "Enabling quadrature tests:")

set(TestExe PolynomialQuadratureTests)

# Add target for all tests
add_executable(${TestExe} QuadratureTests.cpp)

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

include(ConfigureTests)
quicc_add_test(ChebyshevRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(WorlandChebyshevRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(WorlandSphEnergyRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(WorlandCylEnergyRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(WorlandLegendreRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(LegendreRuleTest
  COMMAND ${TestExe}
  KEYWORD Rule
  )
quicc_add_test(JacobiRuleTest
 COMMAND ${TestExe}
 KEYWORD Rule
 ULP 650
 )

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Polynomial
  FILENAME "Quadrature.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)
