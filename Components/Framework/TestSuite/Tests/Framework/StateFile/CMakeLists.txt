message(VERBOSE "Enabling Io/Variable tests:")

set(TestExe FrameworkStateFileTests)

# Add target for all tests
add_executable(${TestExe} StateFileTests.cpp)

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
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Framework
  FILENAME "StateFile.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

include(ConfigureTests)

add_custom_command(TARGET ${TestExe} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E remove state*.hdf5
  COMMAND ${CMAKE_COMMAND} -E copy _refdata/Framework/StateFile/state*.hdf5 .
  WORKING_DIRECTORY ${QUICC_WORK_DIR}
  )

set(_schemes
  WLFl
  WLFm
  SLFl
  SLFm
  TFF
  )
foreach(_scheme ${_schemes})
  foreach(_alg serial)
    foreach(_mode read write)
      quicc_add_test(${_scheme}:${_mode}
        COMMAND ${TestExe}
        KEYWORD Scheme:Mode
        ID 0
        OPTIONS algorithm:${_alg}
      )
    endforeach()
    add_test(
      NAME ${TestExe}_${_scheme}_validate_algorithm-${_alg}
      COMMAND h5diff -q state0000.hdf5 state_${_scheme}.hdf5
      WORKING_DIRECTORY "${QUICC_WORK_DIR}"
      )
  endforeach()
endforeach()

#
# profiling tests
#
set(_algs serial)
if(QUICC_USE_MPI)
  list(APPEND _algs tubular)
endif()
foreach(_scheme ${_schemes})
  foreach(_alg ${_algs})
    foreach(_mode read write)
      quicc_add_test(${_scheme}:${_mode}
        COMMAND ${TestExe}
        KEYWORD Scheme:Mode
        ID 0
        OPTIONS algorithm:${_alg}
        STEPS 10
        PERFONLY
      )
    endforeach()
  endforeach()
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
