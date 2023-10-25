message(VERBOSE "Enabling Timestep/RungeKuttaCB tests:")

set(TestExe TimestepRkcbTests)

# Add target for all tests
add_executable(${TestExe} TimestepRkcbTests.cpp)

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
  ImExRKCB2
  ImExRKCB3b
  ImExRKCB3c
  ImExRKCB3d
  ImExRKCB3e
  ImExRKCB3f
  ImExRKCB4
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Scheme
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
