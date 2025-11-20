#
# Utility to add model libraries
#
# target
#     name/path of the model
# TYPES
#     list of model types
# SOURCES_DIRS
#     addition source directories
# EXTRA_LIBS
#     additional libraries to link to
#
function(quicc_add_model target)
  # parse inputs
  set(multiValueArgs TYPES SOURCE_DIRS MODEL_DIRNAME EXTRA_LIBS)
  cmake_parse_arguments(QAM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_model")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAM_TYPES: ${QAM_TYPES}")
  message(DEBUG "QAM_MODEL_DIRNAME: ${QAM_MODEL_DIRNAME}")

  if(NOT QAM_SOURCE_DIRS)
    set(QAM_SOURCE_DIRS Model/)
  endif()
  message(DEBUG "QAM_SOURCE_DIRS: ${QAM_SOURCE_DIRS}")
  message(DEBUG "QAM_EXTRA_LIBS: ${QAM_EXTRA_LIBS}")

  set(_model_dir "Models/${QAM_MODEL_DIRNAME}")
  string(TOLOWER "quicc_${QAM_MODEL_DIRNAME}" _model_lib)

  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

  # Set library visibility
  set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

  # Setup generic model library
  add_library(${_model_lib} "")
  set_target_properties(${_model_lib} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${_model_lib} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    ${PROJECT_BINARY_DIR}/${_model_dir}/Git
    )
  target_link_libraries(${_model_lib} PUBLIC
    QuICC::Framework
    QuICC::DenseSM
    )
  foreach(_extra_lib ${QAM_EXTRA_LIBS})
    target_link_libraries(${_model_lib} PUBLIC
      ${_extra_lib}
      )
  endforeach()

  # Update python files
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Python")
    add_custom_target(${_model_lib}_updatepy)
    add_custom_command(TARGET ${_model_lib}_updatepy POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      "${CMAKE_CURRENT_SOURCE_DIR}/Python"
      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python"
      COMMENT "Copying Python files for ${QAM_MODEL_DIRNAME}"
      VERBATIM
      )
    add_dependencies(${_model_lib} ${_model_lib}_updatepy)
    if(QUICC_CURRENT_UPDATEPY_TARGET)
      add_dependencies(${QUICC_CURRENT_UPDATEPY_TARGET} ${_model_lib}_updatepy)
    endif()
    set(QUICC_CURRENT_UPDATEPY_TARGET "${_model_lib}_updatepy" CACHE STRING "Make dependencies across updatepy" FORCE)
  endif()

  # Generate git hash library
  include(gitUtils/AddGitHashLib)
  string(REPLACE "/" "::" _tgt_nmsp "Model/${target}")
  AddGitHashLib(NAMESPACE ${_tgt_nmsp})
  # Link
  target_link_libraries(${_model_lib} PUBLIC
    "${QUICC_NAMESPACE}${_tgt_nmsp}::GitHash"
  )

  # Create model implementation libraries
  foreach(type ${QAM_TYPES})
    string(TOLOWER "${_model_lib}_${type}" modLib)
    message(DEBUG "modLib: ${modLib}")
    add_library(${modLib} "")
    set_target_properties(${modLib} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${modLib} PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    )
    target_link_libraries(${modLib} PUBLIC
      ${_model_lib}
    )

    string(TOUPPER "QUICC_MODEL_${QAM_MODEL_DIRNAME}_${type}_BACKEND" _modBackend)
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

  foreach(src ${QAM_SOURCE_DIRS})
    add_subdirectory(${src})
  endforeach()

  unset(QAM_TYPES)
  unset(modLib)

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
