#
# Create the Config executable
#
function (quicc_create_config_exe ModelId ModelLib)
  quicc_create_all_exe("${ModelId}" "Config")
  quicc_add_exe("${ModelId}" "Config" "${QUICC_EXE_DIR}/WriteConfig.cpp" "${ModelLib}")
endfunction (quicc_create_config_exe)


#
# Create the State executable
#
function (quicc_create_state_exe ModelId ModelLib)
  quicc_create_all_exe("${ModelId}" "State")
  quicc_add_exe("${ModelId}" "State" "${QUICC_EXE_DIR}/GenerateState.cpp" "${ModelLib}")
endfunction (quicc_create_state_exe)


#
# Create the Run executable
#
function (quicc_create_model_exe ModelId ModelLib)
  quicc_create_all_exe("${ModelId}" "Model")
  quicc_add_exe("${ModelId}" "Model" "${QUICC_EXE_DIR}/RunSimulation.cpp" "${ModelLib}")
endfunction (quicc_create_model_exe)


#
# Create the RunStability executable
#
function (quicc_create_stability_exe ModelId ModelLib)
  if(QUICC_HAVE_STABILITY_SOLVER)
    quicc_create_all_exe("${ModelId}" "Stability")
    quicc_add_exe("${ModelId}" "Stability" "${QUICC_TOOLS_DIR}/Stability/RunStability.cpp" "${ModelLib}"
      EXTRALIBS QuICC::Stability)
  endif()
endfunction (quicc_create_stability_exe)


#
# Create the Visu executable
#
function (quicc_create_visu_exe ModelId ModelLib)
  quicc_create_all_exe("${ModelId}" "Visu")
  quicc_add_exe("${ModelId}" "Visu" "${QUICC_EXE_DIR}/VisualizeState.cpp" "${ModelLib}")
endfunction (quicc_create_visu_exe)


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
function (quicc_add_exe ModelId Postfix ExeSrc ModelLib)
  # parse inputs
  set(multiValueArgs EXTRALIBS)
  cmake_parse_arguments(QAE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

  # Create simple model name
  quicc_model_id2name(ModelName ${ModelId})
  quicc_model_id2cpp(CPPModel ${ModelId})

  # Create new name for executable
  set(ExeName ${ModelName}${Postfix})

  # Add executable to target list
  add_executable(${ExeName} ${ExeSrc})

  # Link to model library
  target_link_libraries(${ExeName}
    ${ModelLib}
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
    RUNTIME_OUTPUT_DIRECTORY "Executables/"
    )
  target_compile_definitions(${ExeName} PRIVATE
    "QUICC_RUNSIM_PATH=${ModelId}"
    "QUICC_RUNSIM_CPPMODEL=${CPPModel}"
    )

  # Install
  install(TARGETS ${ExeName})

  # Show message
  message(VERBOSE "added ${ExeName}")

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction ()
