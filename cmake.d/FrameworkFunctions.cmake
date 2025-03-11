#
# Backport of new target_sources that uses relative paths
#
function(quicc_target_sources target)
  if(POLICY CMP0076)
    # New behavior is available, so just forward to it by ensuring
    # that we have the policy set to request the new behavior, but
    # don't change the policy setting for the calling scope
    cmake_policy(PUSH)
    cmake_policy(SET CMP0076 NEW)
    target_sources(${target} ${ARGN})
    cmake_policy(POP)
    return()
  endif()

  # Must be using CMake 3.12 or earlier, so simulate the new behavior
  unset(_srcList)
  get_target_property(_targetSourceDir ${target} SOURCE_DIR)

  foreach(src ${ARGN})
    if(NOT src STREQUAL "PRIVATE" AND
        NOT src STREQUAL "PUBLIC" AND
        NOT src STREQUAL "INTERFACE" AND
        NOT IS_ABSOLUTE "${src}")
      # Relative path to source, prepend relative to where target was defined
      file(RELATIVE_PATH src "${_targetSourceDir}" "${CMAKE_CURRENT_LIST_DIR}/${src}")
    endif()
    list(APPEND _srcList ${src})
  endforeach()
  target_sources(${target} ${_srcList})
endfunction()


#
# target_sources that forces the use of CUDA if enabled
#
function(quicc_target_cuda_sources toggle target)
  get_property(_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(${toggle} AND "CUDA" IN_LIST _languages)
    unset(_cudaList)
    foreach(src ${ARGN})
      if(NOT src STREQUAL "PRIVATE" AND
          NOT src STREQUAL "PUBLIC" AND
          NOT src STREQUAL "INTERFACE" AND
          NOT IS_ABSOLUTE "${src}")
        set(_abssrc "${src}")
        string(REPLACE ".cpp" ".cu" _cusrc "${src}")
        set(_abscusrc "${_cusrc}")
        string(PREPEND _abssrc "${CMAKE_CURRENT_LIST_DIR}/")
        string(PREPEND _abscusrc "${CMAKE_CURRENT_LIST_DIR}/")
        file(CREATE_LINK "${_abssrc}" "${_abscusrc}")
      else()
        set(_cusrc "${src}")
      endif()
      list(APPEND _cudaList ${_cusrc})
    endforeach()

    # Call standard quicc_target_sources
    quicc_target_sources("${target}" ${_cudaList})
  else()
    quicc_target_sources("${target}" ${ARGN})
  endif()
endfunction()
