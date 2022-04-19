function(quicc_add_framework_transform_cartesianchebyshev_test target test_exe)
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

  message(VERBOSE "Adding ${target} from template")
  set(testfile
    "${PROJECT_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Tests/Framework/Transform/CartesianChebyshev/${target}.cpp")
  set(specN "128")
  set(physN "256")
  set(blockSize "13")
  set(idPairs "{0,1},{1,2},{2,3}")
  set(lower "-1.0")
  set(upper "1.0")
  configure_file(
    "${target}.cpp.in"
    "${testfile}"
  )

  quicc_target_sources(${test_exe} PRIVATE
    ${testfile}
  )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_add_framework_transform_cartesianchebyshev_test)

message(VERBOSE "Enabling CartesianChebyshev tests:")

set(TestExe FrameworkTransformCartesianChebyshevTests)

# Add target for all tests
add_executable(${TestExe} CartesianChebyshevTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite"
)

# Create Cartesian Chebyshev tests targets
quicc_add_framework_transform_cartesianchebyshev_test(InitTest ${TestExe})

file(MAKE_DIRECTORY
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/data/Framework/Transform/CartesianChebyshev")

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
