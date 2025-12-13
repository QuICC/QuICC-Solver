message(VERBOSE "Enabling PhysicalOperators tests:")

set(TestExe FrameworkPhysicalOperatorsTests)

add_executable(${TestExe} 
  details/Helper.cpp
  PhysicalOperatorsTests.cpp
  CrossTests.cpp
  DotTests.cpp
  MassFluxAdvectionTests.cpp
  SphericalBuoyancyTests.cpp
  SphericalCoriolisTests.cpp
  SphericalHeatAdvectionTests.cpp
  SphericalSComponentTests.cpp
  SphericalZComponentTests.cpp
  StreamAdvectionTests.cpp
  StreamHeatAdvectionTests.cpp
  VelocityAdvectionTests.cpp
  VelocityHeatAdvectionTests.cpp
)

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

add_test(
  NAME ${TestExe}
  COMMAND ${TestExe}
  WORKING_DIRECTORY
  "${QUICC_WORK_DIR}"
)
