function(quicc_set_default_git_branch name)
  # parse inputs
  set(oneValueArgs BRANCH TYPE)
  cmake_parse_arguments(QSGB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  string(TOUPPER "${QSGB_TYPE}" _type)
  message(DEBUG "QSGB_TYPE: ${QSGB_TYPE}")
  message(DEBUG "QSGB_BRANCH: ${QSGB_BRANCH}")

  set(QUICC_${_type}_GIT_BRANCH_${name}
    ${QSGB_BRANCH} CACHE STRING "${QSGB_TYPE} repository branch.")
  set(_QUICC_${_type}_GIT_BRANCH_${name} "${QSGB_BRANCH}" CACHE INTERNAL "${QSGB_TYPE} default repository branch.")
endfunction()

function(quicc_fetch_repository name)
  # parse inputs
  set(oneValueArgs TYPE)
  cmake_parse_arguments(QFR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  string(TOUPPER "${QFR_TYPE}" _type)

  message(DEBUG "name: ${name}")

  include(gitUtils/GitHubProtocol)
  set(_default_git_repo_base "${QUICC_GITHUB_PREFIX}QuICC/${QFR_TYPE}-")
  message(VERBOSE "Default ${QFR_TYPE} repository: ${_default_git_repo_base}")
  set(QUICC_${_type}_GIT_BRANCH "main" CACHE STRING "${QFR_TYPE} repository branch.")
  mark_as_advanced(QUICC_${_type}_GIT_BRANCH)

  set(_path "${QFR_TYPE}s")

  set(_git_branch "${QUICC_${_type}_GIT_BRANCH}}")
  include(FetchContent)

  include(SubDirList)

  # absolute type path
  set(_quicc_path_abs ${PROJECT_SOURCE_DIR}/${_path})

  # grab all folders so we can check against them
  # in a case insensitive way
  subdirlist(_all_dirs ${_quicc_path_abs})

  # check if folder exists
  string(TOLOWER ${name} _name_lc)
  set(_name_exist "False")
  foreach(_dir ${_all_dirs})
      message(DEBUG "_dir: ${_dir}")
      message(DEBUG "_name_lc: ${_name_lc}")
      string(TOLOWER ${_dir} _dir_lc)
      if(${_dir_lc} STREQUAL ${_name_lc})
          set(_name_exist "True")
          # use existing directory name
          set(_name ${_dir})
          break()
      endif()
  endforeach()

  set(_name_path ${_quicc_path_abs}/${name})
  message(DEBUG "_name_path: ${_name_path}")

  if(_name_exist)
    # if the folder exists with a valid CMakeLists.txt, add it
    if(NOT EXISTS ${_name_path}/CMakeLists.txt)
      message(FATAL_ERROR "${_name_path} exists, but does not contain a valid ${QFR_TYPE}.")
    endif()
    message(VERBOSE "${name} exists")
    if(QUICC_${_type}_GIT_BRANCH_${name})
        # check if the branch/tag is the same as default
        include(gitUtils/GetGitBranchTag)
        quicc_get_branch(_name_branch PATH ${_name_path})
        set(_repo_branch ${_QUICC_${_type}_GIT_BRANCH_${name}})
        list(FIND _name_branch ${_repo_branch} _pos)
        if(${_pos} LESS 0)
            message(WARNING "${name} default (${_repo_branch}) and existing branch (${_name_branch}) do not match")
        endif()
    endif()
    message(VERBOSE "adding ${QFR_TYPE}..")
    add_subdirectory(${_name_path})
  else()
      # otherwise try to fetch it from repo
      message(VERBOSE "${name} does not exist, trying to fetch it..")

      # check for type specific repo
      if(QUICC_${_type}_GIT_REPO_BASE_${name})
          set(_repo_base ${QUICC_${_type}_GIT_REPO_BASE_${name}})
      # check for custom repo
      elseif(QUICC_${_type}_GIT_REPO_BASE)
          set(_repo_base ${QUICC_${_type}_GIT_REPO_BASE})
      else()
          set(_repo_base ${_default_git_repo_base})
      endif()
      # check for type specific branch
      message(DEBUG "QUICC_${_type}_GIT_BRANCH_${name}: ${QUICC_${_type}_GIT_BRANCH_${name}}")
      if(QUICC_${_type}_GIT_BRANCH_${name})
          set(_repo_branch ${QUICC_${_type}_GIT_BRANCH_${name}})
      else()
          set(_repo_branch ${_git_branch})
      endif()
      set(_repo ${_repo_base}${name})
      message(DEBUG "_repo_base: ${_repo_base}")
      message(DEBUG "_repo_branch: ${_repo_branch}")
      message(DEBUG "_repo: ${_repo}")

      if(NOT CMAKE_MESSAGE_LOG_LEVEL STREQUAL "VERBOSE" AND
         NOT CMAKE_MESSAGE_LOG_LEVEL STREQUAL "DEBUG")
          set(_quiet QUIET)
      endif()

      FetchContent_Declare(
        ${name}
        ${_quiet}
        GIT_REPOSITORY ${_repo}
        GIT_TAG ${_repo_branch}
        GIT_PROGRESS TRUE
        SOURCE_DIR ${_name_path}
        SUBBUILD_DIR ${CMAKE_BINARY_DIR}/${_path}/${name}
        BINARY_DIR ${CMAKE_BINARY_DIR}/${_path}/${name}
      )

      FetchContent_MakeAvailable(${name})
  endif()

  mark_as_advanced_all(FETCHCONTENT)
endfunction()
