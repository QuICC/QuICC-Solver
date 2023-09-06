message(VERBOSE "Enabling PyQuICC tests:")

set(TestExe PyQuICCTests)

# Add target for all tests
add_executable(${TestExe} PyQuICCTests.cpp)

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

# Create PyQuICC tests targets
include(ConfigureTests)
set(_tests
  ToolsTest
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD notused
  )
endforeach()
