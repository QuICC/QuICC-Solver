message(VERBOSE "Enabling LoadSplitter tests:")

set(TestExe FrameworkLoadSplitterTests)

# Add target for all tests
add_executable(${TestExe} LoadSplitterTests.cpp)

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
  COMPONENT Framework
  FILENAME "LoadSplitter.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Test configuratio tool
include(ConfigureTests)

# Spatial scheme list
set(_tests
  WLFl
  SLFl
  WLFm
  SLFm
  TFF
  )

# Serial tests
foreach(_test ${_tests})
  foreach(_db 104 106 108)
    foreach(_st 0 1 2)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:serial np:1 db:${_db} stage:${_st}
      )
    endforeach()
  endforeach()
endforeach()

# Single1D tests
foreach(_test ${_tests})
  foreach(_db 106)
    foreach(_st 0 1 2)
      foreach(_alg "single1d")
        foreach(_np 4 16 42)
          quicc_add_test(${_test}
            COMMAND ${TestExe}
            KEYWORD Scheme
            OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st}
          )
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

# Single2D tests
foreach(_test ${_tests})
  foreach(_db 104)
    foreach(_st 0 1 2)
      foreach(_alg "single2d")
        foreach(_np 4 6 8)
          quicc_add_test(${_test}
            COMMAND ${TestExe}
            KEYWORD Scheme
            OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st}
          )
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

# Tubular tests
foreach(_test ${_tests})
  foreach(_db 106)
    foreach(_st 0 1 2)
      foreach(_alg "tubular")
        foreach(_np 6 128 288)
          quicc_add_test(${_test}
            COMMAND ${TestExe}
            KEYWORD Scheme
            OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st}
          )
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

# Tubular tests 2**8
foreach(_test ${_tests})
  foreach(_db 108)
    foreach(_st 0 1 2)
      foreach(_alg "tubular")
        foreach(_np 144)
          quicc_add_test(${_test}
            COMMAND ${TestExe}
            KEYWORD Scheme
            OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st}
          )
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
