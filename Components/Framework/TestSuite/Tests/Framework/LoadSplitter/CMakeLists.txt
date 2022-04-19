message(VERBOSE "Enabling LoadSplitter tests:")

set(TestExe FrameworkLoadSplitterTests)

# Add target for all tests
add_executable(${TestExe} LoadSplitterTests.cpp)

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
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
set(_tests
  Tubular::WLFl
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Scheme
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
