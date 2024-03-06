message(VERBOSE "Enabling Bessel tests:")

set(TestExe SparseSMBesselTests)

# Add target for all tests
add_executable(${TestExe} BesselTests.cpp)

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
  COMPONENT SparseSM
  FILENAME "Bessel.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Create Bessel basis tests targets
include(ConfigureTests)
set(_tests
  SphLapl:Value:sparse
  SphLapl2:Value:sparse
  SphLapl:Insulating:sparse
  SphLapl2:Insulating:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 0 1 2 3 4
    ULP 1800
  )
endforeach()

# Create Bessel basis tests targets
set(_tests
  SphLapl:Value:sparse
  SphLapl2:Value:sparse
  SphLapl:Insulating:sparse
  SphLapl2:Insulating:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 10 11 12 13 14
    ULP 1800
  )
endforeach()

# Create Bessel basis tests targets
set(_tests
  SphLapl:Value:sparse
  SphLapl:Insulating:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 20 21 22 23 24 25 26 27
    ULP 1800
  )
endforeach()

# Create Bessel basis tests targets
set(_tests
  SphLapl2:Value:sparse
  SphLapl2:Insulating:sparse
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:BType:type
    IDS 20 21 22 23 24 25 26
    ULP 1800
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
