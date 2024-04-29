#
# Use multiple precision computation for inititialisation.
#
option(QUICC_MULTPRECISION "Enable multiple precision computations?" OFF)
if(QUICC_MULTPRECISION)
  quicc_create_option(NAME QUICC_MPBACKEND
      OPTS "boost" "gmp" "mpfr" "quad"
      LABEL "Multiple precision backend"
  )

  if(NOT QUICC_MPBACKEND STREQUAL "boost")
    message(FATAL_ERROR "find gmp/mpfr/quad library needs to be implemented")
  endif()

  include(QuICCFindBoostRequired)
  quicc_find_boost_required()
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
    INTERFACE
      Boost::headers)

  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
    INTERFACE
      Boost::headers)

  target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}
    INTERFACE "-DQUICC_MULTPRECISION")
  quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB}
    INTERFACE OPTION QUICC_MPBACKEND)

  if(NOT QUICC_MPBACKEND STREQUAL "quad")
    if(NOT DEFINED QUICC_MULTPRECISION_DIGITS)
      set(QUICC_MULTPRECISION_DIGITS 50)
    endif()
    set(QUICC_MULTPRECISION_DIGITS ${QUICC_MULTPRECISION_DIGITS} CACHE STRING "Multiple precision digits" FORCE)
    target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}
      INTERFACE "-DQUICC_MULTPRECISION_DIGITS=${QUICC_MULTPRECISION_DIGITS}")
  endif()

  message(STATUS "Multiple precision computations are enabled")

endif()
