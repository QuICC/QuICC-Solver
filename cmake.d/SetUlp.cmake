function(quicc_set_ulp var)
  # parse inputs
  set(oneValueArgs )
  set(multiValueArgs ULP MPULP)
  cmake_parse_arguments(QSU "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_set_ulp")
  message(DEBUG "var: ${var}")
  message(DEBUG "QSU_ULP: ${QSU_ULP}")
  if(NOT QSU_MPULP)
    set(QSU_MPULP ${QSU_ULP})
  endif()
  message(DEBUG "QSU_MPULP: ${QSU_MPULP}")

  if(QUICC_MULTPRECISION)
    set(_ulp ${QSU_MPULP})
  else()
    set(_ulp ${QSU_ULP})
  endif()

  set(${var} ${_ulp} PARENT_SCOPE)

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()

