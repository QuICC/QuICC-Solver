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

# Spatial scheme list for uniform truncation
set(_tests
  WLFl
  SLFl
  WLFm
  SLFm
  TFF
  )

# Serial tests
set(_alg "serial")
set(_np 1)
foreach(_test ${_tests})
  foreach(_db 104 106 108)
    foreach(_st 0 1 2 3)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st}
      )
    endforeach()
  endforeach()
endforeach()

# Single1D tests
set(_alg "single1d")
set(_db 106)
foreach(_test "TFF")
  foreach(_st 0 1 2 3)
    foreach(_np 4 16 42)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()
foreach(_test "WLFl" "SLFl")
  foreach(_st 1 2)
    foreach(_np 4 16 42)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()
foreach(_test "WLFm" "SLFm")
  foreach(_st 1 2 3)
    foreach(_np 4 16 42)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Single2D tests
set(_alg "single2d")
set(_db 104)
foreach(_test ${_tests})
  foreach(_st 0 1 3)
    foreach(_np 4 6 8)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Tubular tests
set(_alg "tubular")
set(_db 106)
foreach(_test ${_tests})
  foreach(_st 0 1 2 3)
    foreach(_np 6 128 288)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Tubular tests 2**8
set(_alg "tubular")
set(_db 108)
foreach(_test ${_tests})
  foreach(_st 0 1 2)
    foreach(_np 144)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Full decomposition tests for single1d
set(_alg "single1d")
set(_db 106)
set(_np 14)
foreach(_test ${_tests})
  foreach(_st 0 1 2 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np}
    )
  endforeach()
endforeach()

# Full decomposition tests for single2d
set(_alg "single2d")
set(_db 106)
set(_np 16)
foreach(_test ${_tests})
  foreach(_st 0 1 2 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np}
    )
  endforeach()
endforeach()

# Full decomposition tests for tubular
set(_alg "tubular")
set(_db 106)
set(_np 36)
foreach(_test ${_tests})
  foreach(_st 0 1 2 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np}
    )
  endforeach()
endforeach()

# Spatial scheme list for triangular truncation
set(_tests
  WLFl
  WLFm
  )

# Serial tests
set(_alg "serial")
set(_np 1)
foreach(_test ${_tests})
  foreach(_db 104 106 108)
    foreach(_st 0 3)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} truncation:triangular
      )
    endforeach()
  endforeach()
endforeach()

# Single1D tests
set(_alg "single1d")
set(_db 106)
foreach(_test "WLFm")
  foreach(_st 3)
    foreach(_np 4 16 42)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} truncation:triangular checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Single2D tests
set(_alg "single2d")
set(_db 104)
foreach(_test ${_tests})
  foreach(_st 0 3)
    foreach(_np 4 6 8)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} truncation:triangular checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Tubular tests
set(_alg "tubular")
set(_db 106)
foreach(_test ${_tests})
  foreach(_st 0 3)
    foreach(_np 6 128 288)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} truncation:triangular checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Tubular tests 2**8
set(_alg "tubular")
set(_db 108)
foreach(_test ${_tests})
  foreach(_st 0)
    foreach(_np 144)
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} truncation:triangular checkRanks
      )
    endforeach()
  endforeach()
endforeach()

# Full decomposition tests for single1d
set(_alg "single1d")
set(_db 106)
set(_np 14)
foreach(_test ${_tests})
  foreach(_st 0 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np} truncation:triangular
    )
  endforeach()
endforeach()

# Full decomposition tests for single2d
set(_alg "single2d")
set(_db 106)
set(_np 16)
foreach(_test ${_tests})
  foreach(_st 0 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np} truncation:triangular
    )
  endforeach()
endforeach()

# Full decomposition tests for tubular
set(_alg "tubular")
set(_db 106)
set(_np 36)
foreach(_test ${_tests})
  foreach(_st 0 3)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      OPTIONS algorithm:${_alg} np:${_np} db:${_db} stage:${_st} id:${_np} truncation:triangular
    )
  endforeach()
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
