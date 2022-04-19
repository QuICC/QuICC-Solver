message(VERBOSE "Enabling TransformTree tests:")

set(TestExe FrameworkTransformTransformTreeTests)

# Add target for all tests
add_executable(${TestExe} TransformTreeTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite"
)

quicc_target_sources(${TestExe} PRIVATE
  SimpleTreeTest.cpp
  TorPolTreeTest.cpp
  ScalarTreeTest.cpp
)

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
