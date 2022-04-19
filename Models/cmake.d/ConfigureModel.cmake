#
# Utility to add models
#
# target
#     name/path of the model
# TYPES
#     list of model types
# BROKENMPI
#     list of broken MPIALGO for implicit solver
#
function(quicc_add_model target)
  # parse inputs
  set(multiValueArgs TYPES BROKENMPI)
  cmake_parse_arguments(QAM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_model")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAM_TYPES: ${QAM_TYPES}")
  message(DEBUG "QAM_BROKENMPI: ${QAM_BROKENMPI}")

  # Set Model name and library name
  get_filename_component(modName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  set(QUICC_CURRENT_MODEL_DIR "Models/${modName}")
  string(TOLOWER "quicc_${modName}" QUICC_CURRENT_MODEL_LIB)

  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

  # Disable broken mpi backend
  list(FIND QAM_BROKENMPI ${QUICC_MPIALGO} _pos)
  if(NOT _pos EQUAL -1)
    list(REMOVE_ITEM QAM_TYPES "Implicit")
  endif()

  # Set library visibility
  set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

  # Setup generic model library
  add_library(${QUICC_CURRENT_MODEL_LIB} "")
  set_target_properties(${QUICC_CURRENT_MODEL_LIB} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${QUICC_CURRENT_MODEL_LIB} PUBLIC
    include/
    )
  target_link_libraries(${QUICC_CURRENT_MODEL_LIB} PUBLIC
    quicc_framework
    )

  # Create model implementations libraries
  foreach(type ${QAM_TYPES})
    string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${type}" modLib)
    message(DEBUG "modLib: ${modLib}")
    add_library(${modLib} "")
    set_target_properties(${modLib} PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(${modLib} PUBLIC
      include/
      )
    target_link_libraries(${modLib} PUBLIC
      ${QUICC_CURRENT_MODEL_LIB}
      )
  endforeach()

  # Update python files
  add_custom_target(${QUICC_CURRENT_MODEL_LIB}_updatepy)
  add_custom_command(TARGET ${QUICC_CURRENT_MODEL_LIB}_updatepy POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/Python"
    "${CMAKE_BINARY_DIR}/Python"
    COMMENT "Copying Python files for ${modName}"
    VERBATIM
    )
  add_dependencies(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CURRENT_MODEL_LIB}_updatepy)

  add_subdirectory(src)

  # Create model executables
  foreach(type ${QAM_TYPES})
    set(ModelId "${target}/${type}")
    message(DEBUG "ModelId: ${ModelId}")
    string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${type}" modLib)
    quicc_create_config_exe(${ModelId} "${modLib}")
    quicc_create_state_exe(${ModelId} "${modLib}")
    quicc_create_model_exe(${ModelId} "${modLib}")
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
