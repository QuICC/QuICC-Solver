#
# Utility to suppress output of imported projects
#

function(message mode)
    if(NOT QUICC_MESSAGE_QUIET OR
      CMAKE_MESSAGE_LOG_LEVEL STREQUAL "VERBOSE" OR
      CMAKE_MESSAGE_LOG_LEVEL STREQUAL "DEBUG")
        _message(${mode} "${ARGN}")
    endif()
endfunction()
