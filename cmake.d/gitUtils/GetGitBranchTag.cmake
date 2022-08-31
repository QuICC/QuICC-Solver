#
# Utility to get current git branch/tag
#

function(quicc_get_branch branch)
  # parse inputs
  set(oneValueArgs PATH)
  cmake_parse_arguments(QGB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_get_branch")
  message(DEBUG "QGB_PATH: ${QGB_PATH}")
  if(NOT QGB_PATH)
    set(QGB_PATH ${CMAKE_CURRENT_SOURCE_DIR}})
  endif()

  execute_process(
      COMMAND git show -s --pretty=%d HEAD
      WORKING_DIRECTORY ${QGB_PATH}
      OUTPUT_VARIABLE _branch
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message(DEBUG "_branch: ${_branch}")
  # make a list out of the tags
  if( _branch MATCHES "\\((.*)\\)")
      message(DEBUG "${CMAKE_MATCH_1}")
      string(REPLACE " " "" _branch "${CMAKE_MATCH_1}")
      string(REPLACE "," ";" _branch "${_branch}")
  endif()
  message(DEBUG "_branch: ${_branch}")
  # cleanup list
  foreach(_tag IN LISTS _branch)
    message(DEBUG "_tag: ${_tag}")
    # strip tag:
    if( _tag MATCHES "tag:(.*)")
      message(DEBUG "${CMAKE_MATCH_1}")
      set(_tag "${CMAKE_MATCH_1}")
    endif()
    # strip remote
    if( _tag MATCHES "\\/(.*)")
      message(DEBUG "${CMAKE_MATCH_1}")
      set(_tag "${CMAKE_MATCH_1}")
    endif()

    list(APPEND _tag_list ${_tag})
  endforeach()

  message(DEBUG "_tag_list: ${_tag_list}")
  set(${branch} "${_tag_list}" PARENT_SCOPE)

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
