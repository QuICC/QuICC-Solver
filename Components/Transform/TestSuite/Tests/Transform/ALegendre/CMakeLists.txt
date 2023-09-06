message(VERBOSE "Enabling ALegendre tests:")

set(TestExe TransformALegendreTests)

# Add target for all tests
add_executable(${TestExe} ALegendreTests.cpp)

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

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Transform
  FILENAME "ALegendre.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Transform
  FILENAME "ALegendre_distributed.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

add_subdirectory(Projector)
add_subdirectory(Integrator)
add_subdirectory(BFLoop)

if(QUICC_USE_KOKKOS)
  target_link_libraries(${TestExe} QuICC::Kokkos)
  quicc_fix_target_sources_kokkos(${TestExe})
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
