function(quicc_add_framework_transform_sphereworland_test target test_exe)
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

  message(VERBOSE "Adding ${target} from template")
  set(testfile
    "${PROJECT_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Tests/Framework/Transform/SphereWorland/${target}.cpp")
  set(specN "128")
  set(physN "256")
    set(l "1,5,20,128")
  configure_file(
    "${target}.cpp.in"
    "${testfile}"
  )

  quicc_target_sources(${test_exe} PRIVATE
    ${testfile}
  )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_add_framework_transform_sphereworland_test)

message(VERBOSE "Enabling SphereWorland tests:")

set(TestExe FrameworkTransformSphereWorlandTests)

# Add target for all tests
add_executable(${TestExe} SphereWorlandTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite"
)

# Create Worland tests targets
quicc_add_framework_transform_sphereworland_test(InitTest ${TestExe})

file(MAKE_DIRECTORY
  "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/data/Framework/Transform/SphereWorland")

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
