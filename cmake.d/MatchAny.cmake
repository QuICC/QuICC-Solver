#
# Utility function to check whether any CMake variable that matches NAME 
# contains a given STRING exactly
#
function(match_any)
    # parse inputs
    set(oneValueArgs NAME STRING)
    cmake_parse_arguments(MFUN "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})
    # get all CMake variables
    get_cmake_property(_variableNames VARIABLES)
    # search for matches
    message(DEBUG "Searching ${MFUN_STRING} in ${MFUN_NAME}*")
    list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
    foreach (_variableName ${_variableNames})
        unset(MATCHED)
        message(DEBUG "_variableName: ${_variableName}: ${${_variableName}}")
        # if empty var skip
        list(LENGTH ${_variableName} _len)
        message(DEBUG "_len: ${_len}")
        if(${_len} LESS 1)
            continue()
        endif()
        # check name
        message(DEBUG "_variableName: ${_variableName}")
        string(REGEX MATCH ${MFUN_NAME} MATCHED ${_variableName})
        if (NOT MATCHED)
            continue()
        endif()
        # matched
        message(DEBUG "Match: ${_variableName}")
        list(FIND ${_variableName} "${MFUN_STRING}" _pos)
        message(DEBUG "_pos: ${_pos}")
        if(${_pos} GREATER_EQUAL 0)
            message(DEBUG "Positive Match1: ${_variableName} == ${${_variableName}}" )
            set(${MFUN_STRING}_IS_USED "True" PARENT_SCOPE)
        endif()
    endforeach()
    list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
