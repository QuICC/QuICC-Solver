#
# Utilities to get and cache git hash
#
message(DEBUG "SetGitHash")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

function(GitWriteHash git_hash COMPONENT)
    set(_file "${CMAKE_CURRENT_BINARY_DIR}/${HASH_FILE}")
    string(TOUPPER "${HASH_FILE}" HFUP)
    string(REPLACE "." "_" HFUP ${HFUP})
    string(REPLACE "/" "_" HFUP ${HFUP})
    string(PREPEND HFUP "QUICC_")
    message(DEBUG "HFUP: ${HFUP}")
    message(DEBUG "write _file: ${_file} with ${git_hash}")
    string(REPLACE "/" ";" _list "${COMPONENT}")
    # Generate namespace opening and closing statements
    foreach(_comp IN LISTS _list)
      message(DEBUG ${_comp})
        STRING(APPEND _header_ns_open
"namespace ${_comp}\n\
{\n\
")
        STRING(APPEND _header_ns_close
          "} // namespace ${_comp}\n\
")
    endforeach()
    # Top part of file
    set(_header_file
    "#ifndef ${HFUP}\n\
#define ${HFUP}\n\n\
namespace QuICC\n\
{\n\
")
    # Opening namespaces
    string(APPEND _header_file "${_header_ns_open}")
    string(APPEND _header_file
"    const char gitHash[10] = \"${git_hash}\"\;\n\
")
    # Closing namespaces
    string(APPEND _header_file "${_header_ns_close}")
    # Bottom part of file
    string(APPEND _header_file
"} // namespace QuICC\n\
#endif // ${HFUP}\n\
")
    message(DEBUG ${_header_file})
    file(WRITE ${_file} ${_header_file})
endfunction()

function(GitReadHash git_hash)
    set(_file "${CMAKE_CURRENT_BINARY_DIR}/${HASH_FILE}")
    message(DEBUG "read _file: ${_file}")
    if (EXISTS ${_file})
        message(DEBUG "_file exists")
        file(READ ${_file} _file_content)
        string(REGEX MATCH "\"[0-9,a-z]*\"" _header_hash ${_file_content})
        string(REPLACE "\"" "" _header_hash ${_header_hash})
        message(DEBUG "_header_hash: ${_header_hash}")
    else()
        message(DEBUG "_file does not exist")
        set(_header_hash "invalid")
    endif ()
    set(${git_hash} ${_header_hash} PARENT_SCOPE)
endfunction()

function(CheckGitHash)
    # Get the latest abbreviated commit hash of the working branch
    execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    message(DEBUG "GIT_HASH: ${GIT_HASH}")
    message(DEBUG "HASH_FILE: ${HASH_FILE}")
    message(DEBUG "COMPONENT: ${COMPONENT}")

    GitReadHash(GIT_HASH_CACHE)
    message(DEBUG "GIT_HASH_CACHE: ${GIT_HASH_CACHE}")

    if (NOT ${GIT_HASH} STREQUAL ${GIT_HASH_CACHE})
        message(DEBUG "git hash does not match")
        # Overwrite
        GitWriteHash(${GIT_HASH} ${COMPONENT})
    else()
        message(DEBUG "git hash matches")
    endif ()

endfunction()

function(SetGitHash HASH_FILE)
    # parse inputs
    set(oneValueArgs COMPONENT)
    cmake_parse_arguments(SGH "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})

    set(COMPONENT "${SGH_COMPONENT}")
    string(REPLACE "/" "" _tgt "${SGH_COMPONENT}")

    # pass arguments in a way that works for the external process
    add_custom_target(${_tgt}_AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_HASH=1
        -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
        -DHASH_FILE=${HASH_FILE}
        -DCOMPONENT=${COMPONENT}
        -P ${CMAKE_SOURCE_DIR}/cmake.d/gitUtils/SetGitHash.cmake
        )
    CheckGitHash()
endfunction()

# This is used to run this function from an external cmake process.
if (RUN_CHECK_GIT_HASH)
    CheckGitHash()
    message(DEBUG "running external process")
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)
