###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

#
# Use multiple precision computation for inititialisation.
#
option(QUICC_MULTPRECISION "Enable multiple precision computations?" OFF)
if(QUICC_MULTPRECISION)
  find_package(Boost REQUIRED)

  quicc_create_option(NAME QUICC_MPBACKEND
                    OPTS "boost" "gmp" "mpfr" "quad"
                    LABEL "Multiple precision backend")
  quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC OPTION QUICC_MPBACKEND)

  if(NOT QUICC_LINALG STREQUAL "Eigen")
    message(SEND_ERROR "------->>> Can't use multiple precision computations with selected implementation <<<-------")
  else(NOT QUICC_LINALG STREQUAL "Eigen")
    target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC "-DQUICC_MULTPRECISION")
    if(NOT QUICC_MPBACKEND STREQUAL "boost")
      message(FATAL_ERROR "find gmp/mpfr/quad library needs to be implemented")
    endif()
    if(NOT QUICC_MPBACKEND STREQUAL "quad")
      if(NOT DEFINED QUICC_MULTPRECISION_DIGITS)
        set(QUICC_MULTPRECISION_DIGITS 50)
      endif(NOT DEFINED QUICC_MULTPRECISION_DIGITS)
      set(QUICC_MULTPRECISION_DIGITS ${QUICC_MULTPRECISION_DIGITS} CACHE STRING "Multiple precision digits" FORCE)
      target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC "-DQUICC_MULTPRECISION_DIGITS=${QUICC_MULTPRECISION_DIGITS}")
    endif(NOT QUICC_MPBACKEND STREQUAL "quad")
    message(STATUS "Multiple precision computations are enabled")
  endif(NOT QUICC_LINALG STREQUAL "Eigen")

endif(QUICC_MULTPRECISION)