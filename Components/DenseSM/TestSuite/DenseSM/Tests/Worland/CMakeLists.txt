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
#quicc_fetch_test_reference(
#  ${TestExe}
#  COMPONENT DenseSM
#  FILENAME "Worland.tar.gz"
#  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
#  DATADIR ${QUICC_WORK_DIR}
#)
#
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

include(ConfigureTests)
include(SetUlp)

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
quicc_set_ulp(_ulp
  ULP 550 750 500 1500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 0 1 2 3 4 5 6
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, dipolar tests
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR2D1R1F:dense
  DivR2FD1R1:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 550 750 500 1500
  MPULP 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103
    ULPS ${_ulp}
  )
endforeach()

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
quicc_set_ulp(_ulp
  ULP 550 750 500 1500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 200 201 202 203 204 205 206 207 208 209
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, quadrupolar tests
set(_tests
  DivR1F:dense
  DivR1CF:dense
  DivR1FC:dense
  DivR2D1R1F:dense
  DivR2FD1R1:dense
  DivR1D1DivR1FD1R1:dense
  DivR1D1DivR1D1R1F:dense
  DivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 550 750 500 1500
  MPULP 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 300 301 302 303
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt projection dipolar tests
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
    IDS 0 1 2 3 4 5 6 7 8 9 10 11 12 13
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection dipolar tests
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
    IDS 100 101 102 103 104 105 106 107 108 109 110 111 112 113
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt projection quadrupolar tests
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
    IDS 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser projection quadrupolar tests
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
    IDS 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323
    ULPS ${_ulp}
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
