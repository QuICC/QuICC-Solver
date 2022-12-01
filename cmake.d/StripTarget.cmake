#
# Utility function to strip imported target compile options
#
# _main_tgt
#     name of the target to strip
# LANGUAGE
#     optional language specification, defaults to CXX
#
function(strip_target _main_tgt)
    # parse inputs
    set(oneValueArgs LANGUAGE)
    cmake_parse_arguments(STRPTGT "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})


    # search for matches
    message(DEBUG "strip_target ${_main_tgt}")
    list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
    message(DEBUG "STRPTGT_LANGUAGE: ${STRPTGT_LANGUAGE}")

    # default to CXX
    if(NOT STRPTGT_LANGUAGE)
        set(STRPTGT_LANGUAGE CXX)
        message(DEBUG "STRPTGT_LANGUAGE: ${STRPTGT_LANGUAGE}")
    endif()

    if(TARGET ${_main_tgt})
        # strip target
        get_target_property(_main_tgt_INTERFACE_COMPILE_OPTIONS ${_main_tgt} INTERFACE_COMPILE_OPTIONS)
                    message(DEBUG "_main_tgt_INTERFACE_COMPILE_OPTIONS: ${_main_tgt_INTERFACE_COMPILE_OPTIONS}")
        set_target_properties(${_main_tgt}
            PROPERTIES
                INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:${STRPTGT_LANGUAGE}>:>")
        get_target_property(_main_tgt_INTERFACE_COMPILE_OPTIONS ${_main_tgt} INTERFACE_COMPILE_OPTIONS)
        message(DEBUG "_main_tgt_INTERFACE_COMPILE_OPTIONS: ${_main_tgt_INTERFACE_COMPILE_OPTIONS}")

        get_target_property(_main_tgt_INTERFACE_LINK_LIBRARIES ${_main_tgt} INTERFACE_LINK_LIBRARIES)
        # strip sub-targets
        set(_store_link "")
        foreach(_tgt ${_main_tgt_INTERFACE_LINK_LIBRARIES})
            message(DEBUG "_tgt: ${_tgt}")
            if(TARGET ${_tgt})
                strip_target(${_tgt} LANGUAGE ${STRPTGT_LANGUAGE})
                list(APPEND _store_link ${_tgt})
            else()
                message(DEBUG "${_tgt} is not a target")
                # is this a def or a lib? if so, keep it
                if(_tgt MATCHES "^\\-D|^\\-l|^\\-L|\\.so$|\\.a$")
                    message(DEBUG "matched")
                    list(APPEND _store_link ${_tgt})
                endif()
            endif()
        endforeach()
        message(DEBUG "_store_link: ${_store_link}")
        # replace sub targets
        set_target_properties(${_main_tgt}
            PROPERTIES
                INTERFACE_LINK_LIBRARIES "${_store_link}")

    else()
        message(DEBUG "${_main_tgt} is not a target")
    endif()

    list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
