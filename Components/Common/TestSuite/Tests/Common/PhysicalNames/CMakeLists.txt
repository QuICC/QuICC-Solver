function(quicc_add_common_physicalnames_test target test_exe)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp")
    message(STATUS " --- --> Adding ${target} from template")
    set(testfile ${target}.cpp)
  else()
    message(STATUS " --- --> Generating ${target} from template")
    string(REGEX REPLACE "Test" "" PhysicalNames ${target})
    set(testfile
      "${PROJECT_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Tests/Common/PhysicalNames/${PhysicalNames}Test.cpp")
    configure_file(
      "TemplateTest.cpp.in"
      "${testfile}"
    )
  endif()

  quicc_target_sources(${test_exe} PRIVATE
    ${testfile}
  )
endfunction(quicc_add_common_physicalnames_test)

message(STATUS " --> Enabling PhysicalNames tests:")

set(TestExe CommonPhysicalNamesTests)

# Add target for all tests
add_executable(${TestExe} PhysicalNamesTests.cpp)

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
quicc_add_common_physicalnames_test(CodensityTest ${TestExe})
quicc_add_common_physicalnames_test(DensityTest ${TestExe})
quicc_add_common_physicalnames_test(DxMeanTemperatureTest ${TestExe})
quicc_add_common_physicalnames_test(DzMeanTemperatureTest ${TestExe})
quicc_add_common_physicalnames_test(EntropyTest ${TestExe})
quicc_add_common_physicalnames_test(FluctMagneticTest ${TestExe})
quicc_add_common_physicalnames_test(FluctMagneticXTest ${TestExe})
quicc_add_common_physicalnames_test(FluctMagneticYTest ${TestExe})
quicc_add_common_physicalnames_test(FluctMagneticZTest ${TestExe})
quicc_add_common_physicalnames_test(FluctTemperatureTest ${TestExe})
quicc_add_common_physicalnames_test(FluctVelocityTest ${TestExe})
quicc_add_common_physicalnames_test(FluctVelocityXTest ${TestExe})
quicc_add_common_physicalnames_test(FluctVelocityYTest ${TestExe})
quicc_add_common_physicalnames_test(FluctVelocityZTest ${TestExe})
quicc_add_common_physicalnames_test(MagneticTest ${TestExe})
quicc_add_common_physicalnames_test(MagneticXTest ${TestExe})
quicc_add_common_physicalnames_test(MagneticYTest ${TestExe})
quicc_add_common_physicalnames_test(MagneticZTest ${TestExe})
quicc_add_common_physicalnames_test(MeanMagneticTest ${TestExe})
quicc_add_common_physicalnames_test(MeanMagneticXTest ${TestExe})
quicc_add_common_physicalnames_test(MeanMagneticYTest ${TestExe})
quicc_add_common_physicalnames_test(MeanMagneticZTest ${TestExe})
quicc_add_common_physicalnames_test(MeanTemperatureTest ${TestExe})
quicc_add_common_physicalnames_test(MeanVelocityTest ${TestExe})
quicc_add_common_physicalnames_test(MeanVelocityXTest ${TestExe})
quicc_add_common_physicalnames_test(MeanVelocityYTest ${TestExe})
quicc_add_common_physicalnames_test(MeanVelocityZTest ${TestExe})
quicc_add_common_physicalnames_test(PressureTest ${TestExe})
quicc_add_common_physicalnames_test(StreamfunctionTest ${TestExe})
quicc_add_common_physicalnames_test(TemperatureTest ${TestExe})
quicc_add_common_physicalnames_test(VelocityTest ${TestExe})
quicc_add_common_physicalnames_test(VelocityXTest ${TestExe})
quicc_add_common_physicalnames_test(VelocityYTest ${TestExe})
quicc_add_common_physicalnames_test(VelocityZTest ${TestExe})
quicc_add_common_physicalnames_test(VorticityTest ${TestExe})
quicc_add_common_physicalnames_test(VorticityXTest ${TestExe})
quicc_add_common_physicalnames_test(VorticityYTest ${TestExe})
quicc_add_common_physicalnames_test(VorticityZTest ${TestExe})

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
