function(AddGitHashLib)
    # parse inputs
    set(oneValueArgs NAMESPACE)
    cmake_parse_arguments(AGHL "${options}" "${oneValueArgs}"
                          "${multiValueArgs}" ${ARGN})

    message(DEBUG "AGHL_NAMESPACE: ${AGHL_NAMESPACE}")

    # Set path and target
    string(REPLACE "::" "/" _path "${AGHL_NAMESPACE}")
    message(DEBUG "_path: ${_path}")
    string(REPLACE "::" "_" _target "${AGHL_NAMESPACE}_GitHash")
    message(DEBUG "_target: ${_target}")

    set(COMPONENT "${_path}")
    set(HASH_FILE "${_path}/gitHash.hpp")

    include(gitUtils/SetGitHash)
    SetGitHash(${HASH_FILE} COMPONENT ${COMPONENT})

    add_library(${_target} INTERFACE)
    target_include_directories(${_target}
      INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )
      add_dependencies(${_target} ${COMPONENT}_AlwaysCheckGit)

    # Alias
    string(REPLACE "_" "::" _target_al ${_target})
    add_library(${QUICC_NAMESPACE}${_target_al}
      ALIAS ${_target})

    # Export info
    quicc_export_target(${_target}
      COMPONENT ${_target}
      DIRECTORIES
        ${CMAKE_CURRENT_BINARY_DIR}
      FILES_MATCHING_PATTERN "*.hpp"
        )

endfunction()
