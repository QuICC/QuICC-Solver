message(VERBOSE "Enabling ScalarFields tests:")

set(TestExe ScalarFieldsTests)

add_executable(${TestExe} 
  details/Helper.cpp
  ScalarFieldTests.cpp
  FlatScalarFieldTests.cpp
  ViewScalarFieldTests.cpp
)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
  WORKING_DIRECTORY
  "${QUICC_WORK_DIR}"
)
