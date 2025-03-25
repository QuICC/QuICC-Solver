#
# Convert Model ID into Model name
#
function (quicc_model_id2name Name ModelId)
  string(REGEX REPLACE "/" "" ModelName ${ModelId})
  set(${Name} ${ModelName} PARENT_SCOPE)
endfunction (quicc_model_id2name)


#
# Convert Model ID into Model name
#
function (quicc_model_id2cpp Cpp ModelId)
  string(REGEX REPLACE "/" "::" CPPModel ${ModelId})
  set(${Cpp} ${CPPModel} PARENT_SCOPE)
endfunction (quicc_model_id2cpp)


#
# Create target for all main executables
#
function (quicc_create_all_exe ModelId Postfix)
  # Create simple model name
  quicc_model_id2name(ModelName ${ModelId})

  if(NOT TARGET ${ModelName})
    add_custom_target(${ModelName})
    message(VERBOSE "${ModelName}")
  endif()
  add_dependencies(${ModelName} ${ModelName}${Postfix})
endfunction ()


#
# Create executable
#
function (quicc_add_exe ModelId)
  # parse inputs
  set(oneValueArgs POSTFIX EXESOURCE MODELLIB)
  set(multiValueArgs EXTRALIBS)
  cmake_parse_arguments(QAE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QAE_POSTFIX: ${QAE_POSTFIX}")
  message(DEBUG "QAE_EXESOURCE: ${QAE_EXESOURCE}")
  message(DEBUG "QAE_MODELLIB: ${QAE_MODELLIB}")

  # Create simple model name
  quicc_model_id2name(ModelName ${ModelId})
  quicc_model_id2cpp(CPPModel ${ModelId})

  # Create new name for executable
  set(ExeName ${ModelName}${QAE_POSTFIX})

  # Add executable to target list
  add_executable(${ExeName} ${QAE_EXESOURCE})

  # Link to model library
  target_link_libraries(${ExeName}
    ${QAE_MODELLIB}
    )
  # Add extra libraries
  foreach(_lib ${QAE_EXTRALIBS})
    target_link_libraries(${ExeName}
      ${_lib}
      )
  endforeach()

  # Includes
  target_include_directories(${ExeName} PUBLIC
    "include/"
    )

  # Set special properties of target
  set_target_properties(${ExeName} PROPERTIES
    OUTPUT_NAME ${ExeName}
    RUNTIME_OUTPUT_DIRECTORY "./"
    )

  # Install
  install(TARGETS ${ExeName})

  # Show message
  message(VERBOSE "added ${ExeName}")

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction ()

#
# Setup executable
#
function (quicc_setup_exe ModelId)
  # parse inputs
  set(oneValueArgs POSTFIX EXESOURCE MODELLIB)
  set(multiValueArgs EXTRALIBS)
  cmake_parse_arguments(QSE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QSE_POSTFIX: ${QSE_POSTFIX}")
  message(DEBUG "QSE_EXESOURCE: ${QSE_EXESOURCE}")
  message(DEBUG "QSE_MODELLIB: ${QSE_MODELLIB}")

  quicc_create_all_exe("${ModelId}" "${QSE_POSTFIX}")
  quicc_add_exe("${ModelId}"
    POSTFIX "${QSE_POSTFIX}"
    EXESOURCE "${QSE_EXESOURCE}"
    MODELLIB "${QSE_MODELLIB}"
    EXTRALIBS ${QSE_EXTRALIBS})

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction ()
