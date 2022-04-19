function(quicc_add_framework_transform_shellchebyshev_test target test_exe)
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

  message(VERBOSE "Adding ${target} from template")
  set(testfile
    "${PROJECT_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Tests/Framework/Transform/ShellChebyshev/${target}.cpp")
  set(specN "128")
  set(physN "256")
  set(blockSize "13")
  set(idPairs "{0,1},{1,2},{2,3}")
  set(lower "0.5384615384615384")
  set(upper "1.5384615384615385")
  configure_file(
    "${target}.cpp.in"
    "${testfile}"
  )

  quicc_target_sources(${test_exe} PRIVATE
    ${testfile}
  )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_add_framework_transform_shellchebyshev_test)

message(VERBOSE "Enabling ShellChebyshev tests:")

set(TestExe FrameworkTransformShellChebyshevTests)

# Add target for all tests
add_executable(${TestExe} ShellChebyshevTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite"
)

# Create Shell Chebyshev tests targets
quicc_add_framework_transform_shellchebyshev_test(InitTest ${TestExe})

file(MAKE_DIRECTORY
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/data/Framework/Transform/ShellChebyshev")

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
