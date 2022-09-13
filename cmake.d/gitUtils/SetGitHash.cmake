#
# Utilities to get and cache git hash
#
message(DEBUG "SetGitHash")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

function(GitWriteHash git_hash)
    set(_file "${CMAKE_CURRENT_BINARY_DIR}/${HASH_FILE}")
    string(TOUPPER "${HASH_FILE}" HFUP)
    string(REPLACE "." "_" HFUP ${HFUP})
    message(DEBUG "HFUP: ${HFUP}")
    message(DEBUG "write _file: ${_file} with ${git_hash}")
    set(_header_file
    "#ifndef ${HFUP}\n\
#define ${HFUP}\n\n\
namespace QuICC\n\
{\n\
namespace ${SGH_COMPONENT}\n\
{\n\
    const char gitHash[10] = \"${git_hash}\"\;\n\
} // namespace ${SGH_COMPONENT} \n\
} // namespace QuICC \n\
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

    GitReadHash(GIT_HASH_CACHE)
    message(DEBUG "GIT_HASH_CACHE: ${GIT_HASH_CACHE}")

    if (NOT ${GIT_HASH} STREQUAL ${GIT_HASH_CACHE})
        message(DEBUG "git hash does not match")
        # Overwrite
        GitWriteHash(${GIT_HASH})
    else()
        message(DEBUG "git hash matches")
    endif ()

endfunction()

function(SetGitHash HASH_FILE)
    # parse inputs
    set(oneValueArgs COMPONENT)
    cmake_parse_arguments(SGH "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})


    add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_HASH=1
        -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
        -DHASH_FILE=${HASH_FILE}
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
