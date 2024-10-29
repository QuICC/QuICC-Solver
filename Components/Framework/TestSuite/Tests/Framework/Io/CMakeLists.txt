message(VERBOSE "Enabling Framework IO tests:")

set(TestExe FrameworkIoTests)

# Add target for all tests
add_executable(${TestExe} IoTests.cpp)

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

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
#quicc_fetch_test_reference(
#  ${TestExe}
#  COMPONENT Framework
#  FILENAME "Io.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)

include(ConfigureTests)

#add_custom_command(TARGET ${TestExe} POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E remove state*.hdf5
#  COMMAND ${CMAKE_COMMAND} -E copy _refdata/Framework/StateFile/state*.hdf5 .
#  WORKING_DIRECTORY ${QUICC_WORK_DIR}
#  )

quicc_add_test(Dipolarity
  COMMAND ${TestExe}
  KEYWORD Type
  ID 0
)

list(POP_BACK CMAKE_MESSAGE_INDENT)
