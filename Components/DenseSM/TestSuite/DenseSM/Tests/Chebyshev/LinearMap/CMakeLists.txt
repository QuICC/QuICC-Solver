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

# Create Gaunt coefficient, dipolar tests, r^2 prefactor
set(_tests
  RpDivR2D1R1F:dense
  RpDivR2FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
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

# Create Gaunt coefficient, dipolar tests, r^4 prefactor
set(_tests
  RpDivR1D1F:dense
  RpDivR1D1FC:dense
  RpDivR1D1CF:dense
  RpDivR2D1R1F:dense
  RpDivR2D1R1FC:dense
  RpDivR2FD1R1:dense
  RpDivR2CFD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 7 8 9 10 11 12 13
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, dipolar tests, r^2 prefactor
set(_tests
  RpDivR1F:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104 105 106
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, dipolar tests, r^4 prefactor
set(_tests
  RpDivR1F:dense
  RpDivR1FC:dense
  RpDivR1CF:dense
  RpDivR1D1DivR1FD1R1:dense
  RpDivR1D1DivR1D1R1F:dense
  RpDivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 107 108 109 110 111 112 113
    ULPS ${_ulp}
  )
endforeach()

# Create Gaunt coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  RpDivR2D1R1F:dense
  RpDivR2FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500
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

# Create Gaunt coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  RpDivR1D1F:dense
  RpDivR1D1FC:dense
  RpDivR1D1CF:dense
  RpDivR2D1R1F:dense
  RpDivR2D1R1FC:dense
  RpDivR2FD1R1:dense
  RpDivR2CFD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 210 211 212 213 214 215 216 217 218 219
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  RpDivR1F:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 300 301 302 303 304 305 306 307 308 309 310 311
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  RpDivR1F:dense
  RpDivR1FC:dense
  RpDivR1CF:dense
  RpDivR1D1DivR1FD1R1:dense
  RpDivR1D1DivR1D1R1F:dense
  RpDivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 312 313 314 315 316 317 318 319 320 321 322 323
    ULPS ${_ulp}
  )
endforeach()

# Create Quasi-inverse Gaunt coefficient, dipolar tests, r^2 prefactor
set(_tests
  IqRpDivR2D1R1F:dense
  IqRpDivR2FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
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

# Create Quasi-inverse Gaunt coefficient, dipolar tests, r^4 prefactor
set(_tests
  IqRpDivR1D1F:dense
  IqRpDivR1D1FC:dense
  IqRpDivR1D1CF:dense
  IqRpDivR2D1R1F:dense
  IqRpDivR2D1R1FC:dense
  IqRpDivR2FD1R1:dense
  IqRpDivR2CFD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 7 8 9 10 11 12 13
    ULPS ${_ulp}
  )
endforeach()

# Create Quasi-inverse Elsasser coefficient, dipolar tests, r^2 prefactor
set(_tests
  IqRpDivR1F:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 100 101 102 103 104 105 106
    ULPS ${_ulp}
  )
endforeach()

# Create Quasi-inverse Elsasser coefficient, dipolar tests, r^4 prefactor
set(_tests
  IqRpDivR1F:dense
  IqRpDivR1FC:dense
  IqRpDivR1CF:dense
  IqRpDivR1D1DivR1FD1R1:dense
  IqRpDivR1D1DivR1D1R1F:dense
  IqRpDivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 107 108 109 110 111 112 113
    ULPS ${_ulp}
  )
endforeach()

# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  IqRpDivR2D1R1F:dense
  IqRpDivR2FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500
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

# Create Quasi-inverse Gaunt coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  IqRpDivR1D1F:dense
  IqRpDivR1D1FC:dense
  IqRpDivR1D1CF:dense
  IqRpDivR2D1R1F:dense
  IqRpDivR2D1R1FC:dense
  IqRpDivR2FD1R1:dense
  IqRpDivR2CFD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 210 211 212 213 214 215 216 217 218 219
    ULPS ${_ulp}
  )
endforeach()

# Create Elsasser Quasi-inverse coefficient, quadrupolar tests, r^2 prefactor
set(_tests
  IqRpDivR1F:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 300 301 302 303 304 305 306 307 308 309 310 311
    ULPS ${_ulp}
  )
endforeach()

# Create Quasi-inverse Elsasser coefficient, quadrupolar tests, r^4 prefactor
set(_tests
  IqRpDivR1F:dense
  IqRpDivR1FC:dense
  IqRpDivR1CF:dense
  IqRpDivR1D1DivR1FD1R1:dense
  IqRpDivR1D1DivR1D1R1F:dense
  IqRpDivR3D1R1FD1R1:dense
  )
quicc_set_ulp(_ulp
  ULP 500 500 500 500 500 500 500 500 500 500 500 500
  MPULP 25 25 25 25 25 25 25 25 25 25 25 25
        )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Operator:type
    IDS 312 313 314 315 316 317 318 319 320 321 322 323
    ULPS ${_ulp}
  )
endforeach()

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
