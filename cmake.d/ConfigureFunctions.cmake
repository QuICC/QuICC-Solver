#
# Create an option
# If DEFAULT is not provided use first element in OPTS
# If option is provided from cli it has the precendence
# replace quicc_provide_choice and friends
#
function (quicc_create_option)
  # parse inputs
  set(options ADVANCED)
  set(oneValueArgs NAME LABEL DEFAULT)
  set(multiValueArgs OPTS)
  cmake_parse_arguments(QCO "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_create_option")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QCO_NAME: ${QCO_NAME}: ${${QCO_NAME}}")
  message(DEBUG "QCO_OPTS: ${QCO_OPTS}")
  message(DEBUG "QCO_LABEL: ${QCO_LABEL}")
  message(DEBUG "QCO_DEFAULT: ${QCO_DEFAULT}")
  message(DEBUG "QCO_ADVANCED: ${QCO_ADVANCED}")

  # Create string of possible values
  foreach(choice ${QCO_OPTS})
    if(DEFINED _option_labels)
      set(_option_labels "${_option_labels}, ${choice}")
    else(DEFINED _option_labels)
      set(_option_labels "${choice}")
    endif(DEFINED _option_labels)
  endforeach()
  set(_option_labels "Available ${QCO_LABEL}(s): ${_option_labels}")

  # check property of var to honor cli setting
  get_property(_type CACHE ${QCO_NAME} PROPERTY TYPE)
  message(DEBUG "type: ${_type}")
  set(_def_from_cli "FALSE")
  if(DEFINED _type)
    if(${_type} STREQUAL "UNINITIALIZED")
      set(_def_from_cli "TRUE")
    endif()
  endif()
  message(DEBUG "default from cli: ${_def_from_cli}")

  # Create toggle-able options for available entries
  if(NOT DEFINED CACHE{${QCO_NAME}} OR _def_from_cli)
    message(DEBUG "initialize chache var: ${QCO_NAME}")
    # cli
    if(_def_from_cli)
      # check if it is an available options
      list(FIND QCO_OPTS ${${QCO_NAME}} _pos)
      if(_pos EQUAL -1)
        message(WARNING "Provided option ${${QCO_NAME}} for ${QCO_NAME} is not available - ignored")
        message(WARNING "${_option_labels}")
      else()
        set(QCO_DEFAULT "${${QCO_NAME}}")
      endif()
    endif()
    # set default if not provided
    if(NOT QCO_DEFAULT)
      list(GET QCO_OPTS 0 QCO_DEFAULT)
    endif()
    set(${QCO_NAME} ${QCO_DEFAULT} CACHE STRING ${_option_labels} FORCE)
    set_property(CACHE ${QCO_NAME} PROPERTY STRINGS ${QCO_OPTS})
  endif()

  if(QCO_ADVANCED)
    mark_as_advanced(${QCO_NAME})
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  message(VERBOSE "${QCO_NAME}: ${${QCO_NAME}}")
endfunction(quicc_create_option)



#
# Add definition to flags
#
function (quicc_add_definition base)
  # parse inputs
  set(oneValueArgs FORCE)
  cmake_parse_arguments(QAD "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  message(DEBUG "quicc_add_definition")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QAD_FORCE: ${QAD_FORCE}")
  if(QAD_FORCE)
    string(TOUPPER "${QAD_FORCE}" _VALUE)
    string(TOUPPER "${base}_${_VALUE}" def)
    message(DEBUG "def: ${def}")
    add_compile_definitions("${def}")
  else()
    string(TOUPPER "${${base}}" _BASE)
    message(DEBUG "_BASE: ${_BASE}")
    if(NOT ${_BASE} STREQUAL "NONE")
      string(TOUPPER "${base}_${_BASE}" def)
      message(DEBUG "def: ${def}")
      add_compile_definitions("${def}")
    endif()
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction (quicc_add_definition base)


#
# Add definition to flags for a specific target
#
function (quicc_target_add_definition TGT KIND)
  # parse inputs
  set(oneValueArgs OPTION)
  cmake_parse_arguments(QTAD "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  message(DEBUG "quicc_target_add_definition")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "TGT: ${TGT}")
  message(DEBUG "QTAD_OPTION: ${QTAD_OPTION}: ${${QTAD_OPTION}}")
  string(TOUPPER "${${QTAD_OPTION}}" _BASE)
  message(DEBUG "_BASE: ${_BASE}")
  if(NOT ${_BASE} STREQUAL "NONE")
    string(TOUPPER "${QTAD_OPTION}_${_BASE}" def)
    target_compile_definitions(${TGT} ${KIND} ${def})
  endif()
  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction (quicc_target_add_definition)


#
# Export target
#
function (quicc_export_target TGT)
  # parse inputs
  set(oneValueArgs COMPONENT FILES_MATCHING_PATTERN)
  set(multiValueArgs DIRECTORIES)
  cmake_parse_arguments(QET "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  message(DEBUG "quicc_target_add_definition")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "TGT: ${TGT}")
  message(DEBUG "QET_COMPONENT: ${QET_COMPONENT}")
  message(DEBUG "QET_DIRECTORIES: ${QET_DIRECTORIES}")
  message(DEBUG "QET_FILES_MATCHING_PATTERN: ${QET_FILES_MATCHING_PATTERN}")

  # Export info
  install(TARGETS ${TGT}
    EXPORT ${TGT}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  # Set filter to default if not set
  if(NOT QET_FILES_MATCHING_PATTERN)
    set(QET_FILES_MATCHING_PATTERN "*.h*")
    message(DEBUG "QET_FILES_MATCHING_PATTERN: ${QET_FILES_MATCHING_PATTERN}")
  endif()
  # Install header files
  foreach( _dir ${QET_DIRECTORIES})
    install(DIRECTORY ${_dir}
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      COMPONENT ${QET_COMPONENT}
      FILES_MATCHING PATTERN ${QET_FILES_MATCHING_PATTERN}
      )
  endforeach()

  # Export targets
  install(EXPORT ${TGT}Targets
    FILE ${TGT}Targets.cmake
    NAMESPACE ${QUICC_NAMESPACE}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/quicc
  )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction (quicc_export_target)
