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
