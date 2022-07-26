#
# Utility to extend list funtionality to use regex to find a match
#
# list(FIND_REGEX <list> <regex> <output variable>)
#
function(quicc_list)
    if(NOT ${ARGV0} STREQUAL "FIND_REGEX")
        list(${ARGV})
    else()
        list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

        message(DEBUG "ARGV0: ${ARGV0}: ${${ARGV0}}")
        message(DEBUG "ARGC: ${ARGC}")
        # search for matches
        # MODE LIST REGEX OUT
        foreach (_item ${${ARGV1}})
            unset(MATCHED)
            message(DEBUG "_item: ${_item}")
            string(REGEX MATCH "${ARGV2}" MATCHED ${_item})
            if (NOT MATCHED)
                continue()
            else()
                # matched
                message(DEBUG "Match: ${_item}")
                set(${ARGV3} ${_item} PARENT_SCOPE)
                break()
            endif()
        endforeach()
        list(POP_BACK CMAKE_MESSAGE_INDENT)
    endif()
endfunction()
