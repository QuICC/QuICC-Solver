# Check access to GitHub
execute_process(
    COMMAND ssh -o StrictHostKeyChecking=no -T git@github.com
    RESULT_VARIABLE _git_check
    ERROR_QUIET)
message(DEBUG "_git_check: ${_git_check}")
string(FIND ${_git_check} "Permission denied" _pos)
message(DEBUG "_pos: ${_pos}")
if(${_pos} GREATER_EQUAL 0)
    # no ssh access to GitHub
    set(_github_default "https")
else()
    set(_github_default "ssh")
endif()

quicc_create_option(NAME QUICC_GITHUB_PROTOCOL
                    OPTS "https" "ssh"
                    DEFAULT ${_github_default}
                    LABEL "GitHub access protocol.")

if(QUICC_GITHUB_PROTOCOL STREQUAL "https")
    set(_github_protocol_prefix "https://github.com/")
else()
    set(_github_protocol_prefix "git@github.com:")
endif()

set(QUICC_GITHUB_PREFIX ${_github_protocol_prefix})
message(DEBUG ${QUICC_GITHUB_PREFIX})

mark_as_advanced(QUICC_GITHUB_PREFIX)
