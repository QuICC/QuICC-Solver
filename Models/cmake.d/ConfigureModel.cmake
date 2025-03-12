#
# Utility to add models
#
# target
#     name/path of the model
# TYPES
#     list of model types
#
function(quicc_add_model target)
  # parse inputs
  set(multiValueArgs TYPES SOURCE_DIRS)
  cmake_parse_arguments(QAM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_model")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAM_TYPES: ${QAM_TYPES}")

  if(NOT QAM_SOURCE_DIRS)
    set(QAM_SOURCE_DIRS Model/)
  endif()
  message(DEBUG "QAM_SOURCE_DIRS: ${QAM_SOURCE_DIRS}")

  # Set Model name and library name
  get_filename_component(modName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  set(QUICC_CURRENT_MODEL_DIR "Models/${modName}")
  string(TOLOWER "quicc_${modName}" QUICC_CURRENT_MODEL_LIB)

  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

  # Set library visibility
  set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

  # Setup generic model library
  add_library(${QUICC_CURRENT_MODEL_LIB} "")
  set_target_properties(${QUICC_CURRENT_MODEL_LIB} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${QUICC_CURRENT_MODEL_LIB} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    ${PROJECT_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/Git
    )
  target_link_libraries(${QUICC_CURRENT_MODEL_LIB} PUBLIC
    QuICC::Framework
    )

  # Create model implementations libraries
  foreach(type ${QAM_TYPES})
    string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${type}" modLib)
    message(DEBUG "modLib: ${modLib}")
    add_library(${modLib} "")
    set_target_properties(${modLib} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${modLib} PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
      )
    target_link_libraries(${modLib} PUBLIC
      ${QUICC_CURRENT_MODEL_LIB}
      )

    string(TOUPPER "QUICC_MODEL_${modName}_${type}_BACKEND" _modBackend)
    quicc_create_option(
      NAME ${_modBackend}
      OPTS "CPP" "Python"
      LABEL "Backend used for model definition"
      )
    if(${_modBackend} STREQUAL "CPP")
      quicc_target_add_definition(${modLib}
        PUBLIC OPTION ${_modBackend})
    endif()
  endforeach()

  # Update python files
  add_custom_target(${QUICC_CURRENT_MODEL_LIB}_updatepy)
  add_custom_command(TARGET ${QUICC_CURRENT_MODEL_LIB}_updatepy POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/Python"
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python"
    COMMENT "Copying Python files for ${modName}"
    VERBATIM
    )
  add_dependencies(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CURRENT_MODEL_LIB}_updatepy)

  foreach(src ${QAM_SOURCE_DIRS})
    add_subdirectory(${src})
  endforeach()

  # Generate git hash library
  include(gitUtils/AddGitHashLib)
  string(REPLACE "/" "::" _tgt_nmsp "Model/${target}")
  AddGitHashLib(NAMESPACE ${_tgt_nmsp})
  # Link
  target_link_libraries(${QUICC_CURRENT_MODEL_LIB} PUBLIC
    "${QUICC_NAMESPACE}${_tgt_nmsp}::GitHash"
  )


  # Create model executables
  foreach(type ${QAM_TYPES})
    set(ModelId "${target}/${type}")
    message(DEBUG "ModelId: ${ModelId}")
    string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${type}" modLib)
    quicc_create_config_exe(${ModelId} "${modLib}")
    quicc_create_state_exe(${ModelId} "${modLib}")
    quicc_create_model_exe(${ModelId} "${modLib}")
    quicc_create_stability_exe(${ModelId} "${modLib}")
    quicc_create_visu_exe(${ModelId} "${modLib}")
  endforeach()

  unset(QUICC_CURRENT_MODEL_LIB)
  unset(modName)
  unset(QAM_TYPES)
  unset(modLib)

  # Teststuite, not all models have a test
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d/setup/TestSuite.cmake")
    include(setup/TestSuite)
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
