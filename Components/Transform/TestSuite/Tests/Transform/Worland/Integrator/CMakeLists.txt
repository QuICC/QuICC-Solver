message(VERBOSE "Integrators:")

# Set different ULP for normal and MP implementation
include(SetUlp)

# Create tests targets
include(ConfigureTests)

#set(_tests
#  Poly::I4DivR1D1R1_I2:integrator
#  Poly::I6CylLaplh_I4D1R1:integrator
#  Poly::I6DivR1D1R1_I4:integrator
#  )

set(_tests
  Poly::P:integrator
  Poly::P_Zero:integrator
  Poly::R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 8 6 50 28 33
  MPULP 3 2 2 2 2
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests
  Poly::DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 32 51 81 62 66 32 42 53 39 42
  MPULP 2 2 2 2 2 3 3 3 3 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I2:integrator
  Poly::I2_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 13 17 170 59 36
  MPULP 1 4 8 4 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests
  Poly::I2DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 40 47 145 73 86 41 206 191 41 51
  MPULP 2 3 14 4 3 1 10 11 2 3
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I4DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 77 99 599 38 91 74 615 453 48 51
  MPULP 2 11 121 4 4 2 127 93 3 6
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I6DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 53 160 1578 76 51 48 607 1258 66 104
  MPULP 2 11 147 6 8 2 165 92 2 10
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::DivR1D1R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 435 56 54 61 46 605 55 56 51 55
  MPULP 6 2 1 1 1 9 2 1 2 1
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I2DivR1D1R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 20 383 1340 1104 1196 21 214 1017 374 925
  MPULP 2 17 64 42 40 1 10 44 31 24
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Poly::I4DivR1D1R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 73 1262 10831 2375 2005 73 699 6655 1747 2191
  MPULP 2 79 464 121 55 2 58 328 49 51
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests
  Poly::P:integrator
  )
quicc_set_ulp(_ulp
  ULP 111 304 1331 2009 3377
  MPULP 2 2 3 4 5
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 204 205 206 207 208
  )
endforeach()

set(_tests
  Poly::I2:integrator
  )
quicc_set_ulp(_ulp
  ULP 68 150 476 1176
  MPULP 7 6 8 7
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 204 205 206 207
  )
endforeach()

######################
# FFT based algorithm
######################

set(_tests
  Fft::P:integrator
  Fft::P_Zero:integrator
  Fft::R1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 3 3 6 4 5
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests
  Fft::DivR1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 5 2 13 10 10 3 4 7 6 5
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Fft::I2:integrator
  Fft::I2_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 17 22 176 65 31
  MPULP 6 6 27 19 6
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4
  )
endforeach()

set(_tests
  Fft::I2DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 16 23 189 48 78 11 213 205 40 59
  MPULP 10 5 28 13 19 7 54 76 9 7
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Fft::I4DivR1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 55 22 732 59 51 46 518 582 44 22
  MPULP 17 16 102 9 19 11 209 280 11 9
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Fft::DivR1D1R1_Zero:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 11 4 3 3 3 10 8 3 4 3
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Fft::I2DivR1D1R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 10 353 1359 1067 1258 11 193 1111 393 884
  MPULP 4 75 278 95 149 2 57 107 163 97
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

set(_tests
  Fft::I4DivR1D1R1_Zero:integrator
  )
quicc_set_ulp(_ulp
  ULP 44 1213 10190 2217 1803 38 658 6820 1591 2111
  MPULP 6 142 550 227 236 4 175 520 43 121
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 0 1 2 3 4 10 11 12 13 14
  )
endforeach()

#-----------------
# Database tests
#-----------------

set(_tests
  Fft::P:integrator
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS 6 6 9 12 21
    IDS 204 205 206 207 208
  )
endforeach()

set(_tests
  Fft::I2:integrator
  )
quicc_set_ulp(_ulp
  ULP 54 130 472 1145
  MPULP 55 20 26 33
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Algorithm::Operator:type
    ULPS ${_ulp}
    IDS 204 205 206 207
  )
endforeach()
