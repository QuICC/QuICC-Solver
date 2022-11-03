function(quicc_add_common_nondimensional_test target test_exe)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp")
    message(STATUS " --- --> Adding ${target} from template")
    set(testfile ${target}.cpp)
  else()
    message(STATUS " --- --> Generating ${target} from template")
    string(REGEX REPLACE "Test" "" NonDimensional ${target})
    set(testfile
      "${PROJECT_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Tests/Common/NonDimensional/${NonDimensional}Test.cpp")
    configure_file(
      "TemplateTest.cpp.in"
      "${testfile}"
    )
  endif()

  quicc_target_sources(${test_exe} PRIVATE
    ${testfile}
  )
endfunction(quicc_add_common_nondimensional_test)

message(STATUS " --> Enabling NonDimensional tests:")

set(TestExe CommonNonDimensionalTests)

# Add target for all tests
add_executable(${TestExe} NonDimensionalTests.cpp)

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
quicc_add_common_nondimensional_test(AlphaTest ${TestExe})
quicc_add_common_nondimensional_test(BetaTest ${TestExe})
quicc_add_common_nondimensional_test(ChandrasekharTest ${TestExe})
quicc_add_common_nondimensional_test(ChiTest ${TestExe})
quicc_add_common_nondimensional_test(DeltaTest ${TestExe})
quicc_add_common_nondimensional_test(EadyTest ${TestExe})
quicc_add_common_nondimensional_test(EkmanTest ${TestExe})
quicc_add_common_nondimensional_test(ElevatorTest ${TestExe})
quicc_add_common_nondimensional_test(ElsasserTest ${TestExe})
quicc_add_common_nondimensional_test(EpsilonTest ${TestExe})
quicc_add_common_nondimensional_test(EtaTest ${TestExe})
quicc_add_common_nondimensional_test(FastMeanTest ${TestExe})
quicc_add_common_nondimensional_test(GammaTest ${TestExe})
quicc_add_common_nondimensional_test(HeatingTest ${TestExe})
quicc_add_common_nondimensional_test(IotaTest ${TestExe})
quicc_add_common_nondimensional_test(KappaTest ${TestExe})
quicc_add_common_nondimensional_test(LambdaTest ${TestExe})
quicc_add_common_nondimensional_test(Lower1dTest ${TestExe})
quicc_add_common_nondimensional_test(Lower2dTest ${TestExe})
quicc_add_common_nondimensional_test(Lower3dTest ${TestExe})
quicc_add_common_nondimensional_test(MagneticEkmanTest ${TestExe})
quicc_add_common_nondimensional_test(MagneticPrandtlTest ${TestExe})
quicc_add_common_nondimensional_test(MagneticReynoldsTest ${TestExe})
quicc_add_common_nondimensional_test(ModifiedElsasserTest ${TestExe})
quicc_add_common_nondimensional_test(MuTest ${TestExe})
quicc_add_common_nondimensional_test(NuTest ${TestExe})
quicc_add_common_nondimensional_test(OmegaTest ${TestExe})
quicc_add_common_nondimensional_test(OmicronTest ${TestExe})
quicc_add_common_nondimensional_test(PhiTest ${TestExe})
quicc_add_common_nondimensional_test(PiTest ${TestExe})
quicc_add_common_nondimensional_test(PoincareTest ${TestExe})
quicc_add_common_nondimensional_test(PrandtlTest ${TestExe})
quicc_add_common_nondimensional_test(PsiTest ${TestExe})
quicc_add_common_nondimensional_test(RayleighTest ${TestExe})
quicc_add_common_nondimensional_test(RescaledTest ${TestExe})
quicc_add_common_nondimensional_test(RhoTest ${TestExe})
quicc_add_common_nondimensional_test(RobertsTest ${TestExe})
quicc_add_common_nondimensional_test(RossbyTest ${TestExe})
quicc_add_common_nondimensional_test(RRatioTest ${TestExe})
quicc_add_common_nondimensional_test(SigmaTest ${TestExe})
quicc_add_common_nondimensional_test(TauTest ${TestExe})
quicc_add_common_nondimensional_test(TaylorTest ${TestExe})
quicc_add_common_nondimensional_test(ThetaTest ${TestExe})
quicc_add_common_nondimensional_test(Upper1dTest ${TestExe})
quicc_add_common_nondimensional_test(Upper2dTest ${TestExe})
quicc_add_common_nondimensional_test(Upper3dTest ${TestExe})
quicc_add_common_nondimensional_test(UpsilonTest ${TestExe})
quicc_add_common_nondimensional_test(XiTest ${TestExe})
quicc_add_common_nondimensional_test(ZetaTest ${TestExe})

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
)
