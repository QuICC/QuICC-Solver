message(VERBOSE "Enabling Worland tests:")

set(TestExe DenseSMWorlandTests)

# Add target for all tests
add_executable(${TestExe} WorlandTests.cpp)

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
  COMPONENT DenseSM
  FILENAME "Worland.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Generate Iq tests
function(run_iq_tests tag)
  # parse inputs
  set(oneValueArgs ULP MPULP)
  set(multiValueArgs TESTS IDS QIS MIN_QIS)
  cmake_parse_arguments(RIT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT RIT_ULP)
    set(RIT_ULP 200)
  endif()
  if(NOT RIT_MPULP)
    set(RIT_MPULP 120)
  endif()
  if(NOT RIT_MIN_QIS)
    foreach(_t ${RIT_TESTS})
      LIST(GET RIT_QIS 0 _min)
      list(APPEND RIT_MIN_QIS ${_min})
    endforeach()
  endif()
  message(DEBUG "run_iq_tests TAG: ${tag}")
  message(DEBUG "RIT_TESTS: ${RIT_TESTS}")
  message(DEBUG "RIT_IDS: ${RIT_IDS}a")
  message(DEBUG "RIT_QIS: ${RIT_QIS}")
  message(DEBUG "RIT_ULP: ${RIT_ULP}")
  message(DEBUG "RIT_MPULP: ${RIT_MPULP}")

  foreach(_qi ${RIT_QIS})
    set(_qids )
    set(_ulps )
    set(_mpulps )
    foreach(_id ${RIT_IDS})
      math(EXPR _qid "${_qi}*1000 + ${_id}")
      list(APPEND _qids ${_qid})
      list(APPEND _ulps ${RIT_ULP})
      list(APPEND _mpulps ${RIT_MPULP})
    endforeach()
    quicc_set_ulp(_ulp
      ULP ${_ulps}
      MPULP ${_mpulps}
            )
    foreach(_test _min IN ZIP_LISTS RIT_TESTS RIT_MIN_QIS)
      if("${_qi}" EQUAL 0 OR "${_qi}" GREATER_EQUAL "${_min}")
        if("${_qi}" GREATER_EQUAL 1)
          set(_testname "Iq${_test}")
        else()
          set(_testname "${_test}")
        endif()
        quicc_add_test(${_testname}
          COMMAND ${TestExe}
          KEYWORD Operator:type
          IDS ${_qids}
          ULPS ${_ulp}
        )
      endif()
    endforeach()
  endforeach()
endfunction()
set(_qis 0)

# Create Gaunt coefficient, dipolar tests
set(_tests
  DivR1D1F:dense
  DivR1D1CF:dense
  DivR1D1FC:dense
  DivR2D1R1F:dense
  DivR2CFD1R1:dense
  DivR2D1R1FC:dense
  DivR2FD1R1:dense
  )
set (_ids 0 1 2 3 4 5 6)
run_iq_tests("Gaunt, dipolar" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 2100)

# Create Elsasser coefficient, dipolar tests
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
set (_min_qis 1 1 2 2 1 1)
set (_ids 100 101 102 103 104 105 106)
run_iq_tests("Elsasser, dipolar" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 1100)

# Create Gaunt coefficient, quadrupolar tests
set(_tests
  DivR1D1F:dense
  DivR1D1CF:dense
  DivR1D1FC:dense
  DivR2D1R1F:dense
  DivR2CFD1R1:dense
  DivR2D1R1FC:dense
  DivR2FD1R1:dense
  )
set (_min_qis 1 1 3 1 1 2 1)
set (_ids 200 201 202 203 204 205 206 207 208 209)
run_iq_tests("Gaunt, quadrupolar" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 2000)

# Create Elsasser coefficient, quadrupolar tests
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
set (_min_qis 1 1 2 2 1 1)
set (_ids 300 301 302 303 304 305 306 307 308 309 310 311)
run_iq_tests("Elsasser, quadrupolar" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 1900)

# Generate Iq proj tests
function(run_iq_proj_tests tag)
  # parse inputs
  set(oneValueArgs ULP MPULP)
  set(multiValueArgs TESTS IDS QIS MIN_QIS)
  cmake_parse_arguments(RIPT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT RIPT_ULP)
    set(RIPT_ULP 1500)
  endif()
  if(NOT RIPT_MPULP)
    set(RIPT_MPULP 220)
  endif()
  if(NOT RIPT_MIN_QIS)
    foreach(_t ${RIPT_TESTS})
      LIST(GET RIPT_QIS 0 _min)
      list(APPEND RIPT_MIN_QIS ${_min})
    endforeach()
  endif()
  message(DEBUG "run_iq_tests TAG: ${tag}")
  message(DEBUG "RIPT_TESTS: ${RIPT_TESTS}")
  message(DEBUG "RIPT_IDS: ${RIPT_IDS}a")
  message(DEBUG "RIPT_QIS: ${RIPT_QIS}")
  message(DEBUG "RIPT_ULP: ${RIPT_ULP}")
  message(DEBUG "RIPT_MPULP: ${RIPT_MPULP}")

  foreach(_qi ${RIPT_QIS})
    set(_qids )
    set(_ulps )
    set(_mpulps )
    foreach(_id ${RIPT_IDS})
      math(EXPR _qid "${_qi}*1000 + ${_id}")
      list(APPEND _qids ${_qid})
      list(APPEND _ulps ${RIPT_ULP})
      list(APPEND _mpulps ${RIPT_MPULP})
    endforeach()
    quicc_set_ulp(_ulp
      ULP ${_ulps}
      MPULP ${_mpulps}
            )
    foreach(_test _min IN ZIP_LISTS RIPT_TESTS RIPT_MIN_QIS)
      if("${_qi}" EQUAL 0 OR "${_qi}" GREATER_EQUAL "${_min}")
        if("${_qi}" GREATER_EQUAL 1)
          set(_testname "${_test}")
        else()
          set(_testname "${_test}")
        endif()
        quicc_add_test(${_testname}
          COMMAND ${TestExe}
          KEYWORD Operator:type
          IDS ${_qids}
          ULPS ${_ulp}
        )
      endif()
    endforeach()
  endforeach()
endfunction()

# Create Gaunt projection dipolar tests
set(_qis 0)
set(_tests
  ProjCurlPolCrossPol:dense
  )
set (_min_qis 1)
set (_ids 0 2 4 6 8 10 12)
run_iq_proj_tests("Gaunt projection, dipolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})
set (_min_qis 1)
set (_ids 1 3 5 7 9 11 13)
run_iq_proj_tests("Gaunt projection, dipolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

set(_tests
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
set (_min_qis 1 1 1)
set (_ids 0 2 4 6 8 10 12)
run_iq_proj_tests("Gaunt projection, dipolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 6500)
set (_min_qis 1 1 3)
set (_ids 1 3 5 7 9 11 13)
run_iq_proj_tests("Gaunt projection, dipolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 6500)

# Create Gaunt projection dipolar tests
set(_qis 2 3)
set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
set (_min_qis 1 1 1 1)
set (_ids 0 2 4 6 8 10 12)
run_iq_proj_tests("Gaunt projection, dipolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})
set (_min_qis 1 1 1 3)
set (_ids 1 3 5 7 9 11 13)
run_iq_proj_tests("Gaunt projection, dipolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

set(_qis 0 2 3)
# Create Elsasser projection dipolar tests
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  ProjCurlCPolCrossPol:dense
  ProjCurlCurlTorCrossTor:dense
  ProjCurlCurlPolCrossPol:dense
  ProjCurlCurlCPolCrossTor:dense
  )
set (_min_qis 1 1 1 1 2 1)
set (_ids 100 102 104 106 108 110 112)
run_iq_proj_tests("Elsasser projection, dipolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})
set (_min_qis 1 1 2 1 2 2)
set (_ids 101 103 105 107 109 111 113)
run_iq_proj_tests("Elsasser projection, dipolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Gaunt projection quadrupolar tests
set(_qis 0)
set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
set (_min_qis 1 1 1 1)
set (_ids 200 202 204 206 208 210 212 214 216 218)
run_iq_proj_tests("Gaunt projection, quadrupolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 1600)
set (_min_qis 1 1 1 3)
set (_ids 203 205 207 209 211 213 215 217 219)
run_iq_proj_tests("Gaunt projection, quadrupolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  )
set (_min_qis 1 1 1)
set (_ids 201)
run_iq_proj_tests("Gaunt projection, quadrupolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

set(_tests
  ProjCurlCurlCPolCrossPol:dense
  )
set (_min_qis 3)
set (_ids 201)
run_iq_proj_tests("Gaunt projection, quadrupolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis} ULP 11000)

# Create Gaunt projection quadrupolar tests
set(_qis 2 3)
set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
set (_min_qis 1 1 1 1)
set (_ids 200 202 204 206 208 210 212 214 216 218)
run_iq_proj_tests("Gaunt projection, quadrupolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})
set (_min_qis 1 1 1 3)
set (_ids 201 203 205 207 209 211 213 215 217 219)
run_iq_proj_tests("Gaunt projection, quadrupolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

set(_qis 0 2 3)
# Create Elsasser projection quadrupolar tests
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  ProjCurlCPolCrossPol:dense
  ProjCurlCurlTorCrossTor:dense
  ProjCurlCurlPolCrossPol:dense
  ProjCurlCurlCPolCrossTor:dense
  )

set (_min_qis 1 1 1 1 2 1)
set (_ids 300 302 304 306 308 310 312 314 316 318 320 322)
run_iq_proj_tests("Elsasser projection, quadrupolar, fA" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})
set (_min_qis 1 1 2 1 2 2)
set (_ids 301 303 305 307 309 311 313 315 317 319 321 323)
run_iq_proj_tests("Elsasser projection, quadrupolar, fB" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

list(POP_BACK CMAKE_MESSAGE_INDENT)
