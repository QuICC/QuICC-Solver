#
# Utility to fetch reference data
#
# target
#     name of the test to attach custom command to
# NAME
#     name of test
# FILENAME
#     filename of the archive
# ARCHIVEDIR
#     location where archive is stored
# DATADIR
#     location where data is extracted to ${DATADIR}/_refdata
# GITTAG
#     Git tag of version
# GITURL
#     URL of the git repository
function(quicc_fetch_reference target)
  # parse inputs
  set(oneValueArgs NAME FILENAME ARCHIVEDIR DATADIR GITTAG GITURL)
  set(multiValueArgs )
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_fetch_reference")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAT_NAME: ${QAT_NAME}")
  message(DEBUG "QAT_FILENAME: ${QAT_FILENAME}")
  message(DEBUG "QAT_ARCHIVEDIR: ${QAT_ARCHIVEDIR}")
  message(DEBUG "QAT_DATADIR: ${QAT_DATADIR}")
  message(DEBUG "QAT_GITTAG: ${QAT_GITTAG}")
  message(DEBUG "QAT_GITURL: ${QAT_GITURL}")

  # Fetch data
  include(ExternalData )
  # set where to store
  set(ExternalData_BINARY_ROOT "${QAT_ARCHIVEDIR}")
  message(DEBUG "ExternalData_BINARY_ROOT: ${ExternalData_BINARY_ROOT}")
  # set where to find the data sha
  # set(ExternalData_SOURCE_ROOT "${QAT_ARCHIVEDIR}")
  # message(WARNING "ExternalData_SOURCE_ROOT: ${ExternalData_SOURCE_ROOT}")

  set(ExternalData_URL_TEMPLATES
    "${QAT_GITURL}/-/raw/${QAT_GITTAG}/ref/${QAT_NAME}/${QAT_FILENAME}" )

  # Create unique reference data target for tar
  get_filename_component(_name ${QAT_FILENAME} NAME_WE)
  set(RefTarget "${QAT_NAME}_ReferenceData_${_name}")
  message(DEBUG "${QAT_NAME}_ReferenceData_${_name}: ${RefTarget}")

  if(NOT TARGET ${RefTarget})
    ExternalData_Expand_Arguments(${RefTarget}
      TarPath
      DATA{${QAT_FILENAME}}
      )
    message(DEBUG "TarPath: ${TarPath}")

    if(CMAKE_MESSAGE_LOG_LEVEL STREQUAL "VERBOSE" OR
      CMAKE_MESSAGE_LOG_LEVEL STREQUAL "DEBUG")
      set(_prog ON)
    else()
      set(_prog OFF)
    endif()
    ExternalData_Add_Target(${RefTarget} SHOW_PROGRESS ${_prog})

    set(_refdata "${QAT_DATADIR}/_refdata")
    add_custom_command(TARGET ${RefTarget} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${_refdata}"
      COMMAND tar -xf "${TarPath}" -C "${_refdata}"
      COMMAND find "${_refdata}" -type d | xargs -I{} realpath --relative-to "${_refdata}" {} | xargs -I{} mkdir -p "${QAT_DATADIR}/_data/{}"
      )
  endif()

  # target depends on tar as well
  add_dependencies(${target} ${RefTarget})

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
