#
# Register Models
#
#
# A model can be manually added as a folder under QuICC/Models and added
# to the list QUICC_MODEL
#
# Otherwise the models added to the list QUICC_MODEL are looked up in the model repository.
#
#
#
# The following variables affect the behaviour:
#
# QUICC_MODEL
#    List of models to be built (default: )
#    Note that the list is treated as case insensitive (due to GitHub limitation)
#    The case-ness of the first model registration is retained for the path
#
# QUICC_MODEL_GIT_REPO_BASE
#    Git repository to be used instead of default.
#    default: ${QUICC_GITHUB_PREFIX}QuICC/Model-
#    Note that the full repository must be called `<base repo><model>`
#
# QUICC_MODEL_GIT_BRANCH
#    Git repository branch to be used (default: main)
#
# QUICC_MODEL_GIT_REPO_BASE_<model>
#    Git repository to be used for a specific model (default: QUICC_MODEL_GIT_REPO_BASE for unknown models)
#
# QUICC_MODEL_GIT_BRANCH_<model>
#    Git repository branch to be used for a specific model (default: QUICC_MODEL_GIT_BRANCH)
#

# Default branch for known models
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereDynamo "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqShellDynamo "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereRTC "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqShellRTC "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereTC "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqShellTC "feature/transform_setup")
set(_QUICC_MODEL_GIT_BRANCH_BoussinesqPlaneRBC "main")

# Cache value, may be overwritten manually
set(QUICC_MODEL_GIT_BRANCH_BoussinesqSphereDynamo
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereDynamo} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqShellDynamo
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqShellDynamo} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqSphereRTC
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereRTC} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqShellRTC
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqShellRTC} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqSphereTC
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqSphereTC} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqShellTC
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqShellTC} CACHE STRING "Model repository branch.")
set(QUICC_MODEL_GIT_BRANCH_BoussinesqPlaneRBC
    ${_QUICC_MODEL_GIT_BRANCH_BoussinesqPlaneRBC} CACHE STRING "Model repository branch.")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

# Teststuite
include(setup/ModelTestSuite)

include(gitUtils/GitHubProtocol)
set(_default_model_git_repo_base "${QUICC_GITHUB_PREFIX}QuICC/Model-")
message(VERBOSE "Default model repository: ${_default_model_git_repo_base}")
set(QUICC_MODEL_GIT_BRANCH "main" CACHE STRING "Model repository branch.")
mark_as_advanced(QUICC_MODEL_GIT_BRANCH)

set(QUICC_MODEL_PATH "Models")

include(FetchContent)

set(QUICC_MODEL CACHE STRING "List of models to build.")
message(DEBUG "QUICC_MODEL: ${QUICC_MODEL}")

include(SubDirList)

# absolute models path
set(_quicc_model_path_abs ${PROJECT_SOURCE_DIR}/${QUICC_MODEL_PATH})

# grab all folders so we can check against them
# in a case insensitive way
subdirlist(_all_model_dirs ${_quicc_model_path_abs})

foreach(_model ${QUICC_MODEL})
    message(DEBUG "_model: ${_model}")

    # check if folder exists
    string(TOLOWER ${_model} _model_lc)
    set(_model_exist "False")
    foreach(_dir ${_all_model_dirs})
        message(DEBUG "_dir: ${_dir}")
        message(DEBUG "_model_lc: ${_model_lc}")
        string(TOLOWER ${_dir} _dir_lc)
        if(${_dir_lc} STREQUAL ${_model_lc})
            set(_model_exist "True")
            # use existing model name
            set(_model ${_dir})
            break()
        endif()
    endforeach()

    set(_model_path ${_quicc_model_path_abs}/${_model})
    message(DEBUG "_model_path: ${_model_path}")

    if(_model_exist)
        # if the model folder exists with a valid CMakeLists.txt, add it
        if(NOT EXISTS ${_model_path}/CMakeLists.txt)
            message(FATAL_ERROR "${_model_path} exists, but does not contain a valid model.")
        endif()
        message(VERBOSE "${_model} exists")
        if(QUICC_MODEL_GIT_BRANCH_${_model})
            # check if the branch/tag is the same as default
            include(gitUtils/GetGitBranchTag)
            quicc_get_branch(_model_branch PATH ${_model_path})
            set(_repo_branch ${_QUICC_MODEL_GIT_BRANCH_${_model}})
            list(FIND _model_branch ${_repo_branch} _pos)
            if(${_pos} LESS 0)
                message(WARNING "${_model} default (${_repo_branch}) and existing model branch (${_model_branch}) do not match")
            endif()
        endif()
        message(VERBOSE "adding model..")
        add_subdirectory(${_model_path})

    else()
        # otherwise try to fetch it from model repo
        message(VERBOSE "${_model} does not exist, trying to fetch it..")

        # check for model specific repo
        if(QUICC_MODEL_GIT_REPO_BASE_${_model})
            set(_repo_base ${QUICC_MODEL_GIT_REPO_BASE_${_model}})
        # check for custom repo
        elseif(QUICC_MODEL_GIT_REPO_BASE)
            set(_repo_base ${QUICC_MODEL_GIT_REPO_BASE})
        else()
            set(_repo_base ${_default_model_git_repo_base})
        endif()
        # check for model specific branch
        message(DEBUG "QUICC_MODEL_GIT_BRANCH_${_model}: ${QUICC_MODEL_GIT_BRANCH_${_model}}")
        if(QUICC_MODEL_GIT_BRANCH_${_model})
            set(_repo_branch ${QUICC_MODEL_GIT_BRANCH_${_model}})
        else()
            set(_repo_branch ${QUICC_MODEL_GIT_BRANCH})
        endif()
        set(_repo ${_repo_base}${_model})
        message(DEBUG "_repo_base: ${_repo_base}")
        message(DEBUG "_repo_branch: ${_repo_branch}")
        message(DEBUG "_repo: ${_repo}")

        if(NOT CMAKE_MESSAGE_LOG_LEVEL STREQUAL "VERBOSE" AND
           NOT CMAKE_MESSAGE_LOG_LEVEL STREQUAL "DEBUG")
            set(_quiet QUIET)
        endif()

        FetchContent_Declare(
            ${_model}
            ${_quiet}
            GIT_REPOSITORY ${_repo}
            GIT_TAG ${_repo_branch}
            GIT_PROGRESS TRUE
            SOURCE_DIR ${_model_path}
            SUBBUILD_DIR ${CMAKE_BINARY_DIR}/${QUICC_MODEL_PATH}/${_model}
            BINARY_DIR ${CMAKE_BINARY_DIR}/${QUICC_MODEL_PATH}/${_model}
        )

        FetchContent_MakeAvailable(${_model})

    endif()
endforeach()

mark_as_advanced_all(FETCHCONTENT)
