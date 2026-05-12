message(VERBOSE "Enabling Worland tests:")

set(TestExe DenseSMWorlandStencilTests)

# Add target for all tests
add_executable(${TestExe} StencilTests.cpp)

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
#include(FetchTestReference)
#quicc_fetch_test_reference(
#  ${TestExe}
#  COMPONENT DenseSM
#  FILENAME "Worland.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Stencil tests
set(_tests
  OrthogonalValue:dense
  OrthogonalValueD1:dense
  OrthogonalValueD2:dense
  )
foreach(_t ${_tests})
  quicc_add_test(${_t}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0
    ULPS 10
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
