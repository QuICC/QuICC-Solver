message(VERBOSE "Enabling Chebyshev tests:")

set(TestExe DenseSMChebyshevTests)

# Add target for all tests
add_executable(${TestExe} ChebyshevTests.cpp)

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
  FILENAME "ChebyshevLinearMap.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

# Generate IqRp tests
function(run_iqrp_tests tag)
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
  message(DEBUG "run_iqrp_tests TAG: ${tag}")
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
          set(_testname "IqRp${_test}")
        else()
          set(_testname "Rp${_test}")
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



### Tests for the new operators written for the anelastic model
# Create anelastic operators coefficients, dipolar tests, r^2 prefactor

set(_tests
  FD1:dense
  FD2:dense
  FD3:dense
  FD4:dense
  )
set(_ids 0 1 2 3 4 5 6)
set(_qis 0) # the Iq operators are not implemented yet
run_iqrp_tests("Anelastic, dipolar, r^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})

# Create anelastic operators coefficients, dipolar tests, r^4 prefactor
set(_tests
  FD1:dense
  FD2:dense
  FD3:dense
  FD4:dense
  )
set (_min_qis 1 1 1 1) 
set (_ids 7 8 9 10 11 12 13)
run_iqrp_tests("Anelastic, dipolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  FD1:dense
  FD2:dense
  FD3:dense
  FD4:dense
  )
set(_ids 200 201 202 203 204 205 206 207 208 209)
run_iqrp_tests("Anelastic, quadrupolar R^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})


# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  FD1:dense
  FD2:dense
  FD3:dense
  FD4:dense
  )
set (_min_qis 1 1 1 1)
set(_ids 210 211 212 213 214 215 216 217 218 219)
run_iqrp_tests("Gaunt, quadrupolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

### end of tests for the new operators written for the anelastic model



set(_qis 0 1 2 3 4)

# Create Quasi-inverse Gaunt coefficient, dipolar tests, r^2 prefactor
set(_tests
  DivR2D1R1F:dense
  DivR2FD1R1:dense
  )
set(_ids 0 1 2 3 4 5 6)
run_iqrp_tests("Gaunt, dipolar, r^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})

# Create Quasi-inverse Gaunt coefficient, dipolar tests, r^4 prefactor
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
set (_ids 7 8 9 10 11 12 13)
run_iqrp_tests("Gaunt, dipolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Quasi-inverse Elsasser coefficient, dipolar tests, r^2 prefactor
set(_tests
  DivR1F:dense
  )
set(_ids 100 101 102 103 104 105 106)
run_iqrp_tests("Elsasser, dipolar R^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})

# Create Quasi-inverse Elsasser coefficient, dipolar tests, r^4 prefactor
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
set (_min_qis 1 1 2 2 1 1)
set(_ids 107 108 109 110 111 112 113)
run_iqrp_tests("Elsasser, dipolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  DivR2D1R1F:dense
  DivR2FD1R1:dense
  )
set(_ids 200 201 202 203 204 205 206 207 208 209)
run_iqrp_tests("Gaunt, quadrupolar R^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})

# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^4 prefactor
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
set(_ids 210 211 212 213 214 215 216 217 218 219)
run_iqrp_tests("Gaunt, quadrupolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Elsasser Quasi-inverse coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  DivR1F:dense
  )
set(_ids 300 301 302 303 304 305 306 307 308 309 310 311)
run_iqrp_tests("Elsasser, quadrupolar R^2" TESTS ${_tests} IDS ${_ids} QIS ${_qis})

# Create Quasi-inverse Elsasser coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
set (_min_qis 1 1 2 2 1 1)
set(_ids 312 313 314 315 316 317 318 319 320 321 322 323)
run_iqrp_tests("Elsasser, quadrupolar R^4" TESTS ${_tests} IDS ${_ids} QIS ${_qis} MIN_QIS ${_min_qis})

# Create Gaunt projection dipolar tests, r^2 prefactor
set(_tests
  ProjCurlPolCrossPol:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4 5 6 7 8 9 10 11 12 13
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt projection dipolar tests, r^4 prefactor
set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection dipolar tests, r^2
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104 105 106 107 108 109 110 111 112 113
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection dipolar tests, r^4
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  ProjCurlCPolCrossPol:dense
  ProjCurlCurlTorCrossTor:dense
  ProjCurlCurlPolCrossPol:dense
  ProjCurlCurlCPolCrossTor:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt projection quadrupolar tests, r^2 prefactor
set(_tests
  ProjCurlPolCrossPol:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt projection quadrupolar tests, r^4 prefactor
set(_tests
  ProjCurlPolCrossPol:dense
  ProjCurlCurlTorCrossPol:dense
  ProjCurlCurlPolCrossTor:dense
  ProjCurlCurlCPolCrossPol:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection quadrupolar tests, r^2 prefactor
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection quadrupolar tests, r^4 prefactor
set(_tests
  ProjCurlTorCrossPol:dense
  ProjCurlPolCrossTor:dense
  ProjCurlCPolCrossPol:dense
  ProjCurlCurlTorCrossTor:dense
  ProjCurlCurlPolCrossPol:dense
  ProjCurlCurlCPolCrossTor:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
