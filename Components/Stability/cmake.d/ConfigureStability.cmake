#
# Utility to add stability solver for a model
#
# target
#     name/path of the model
# TYPES
#     list of model types
#
include(ModelFunctions)

function(quicc_add_stability target)
  # parse inputs
  set(multiValueArgs TYPES MODEL_DIRNAME)
  cmake_parse_arguments(QAS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_stability")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAS_TYPES: ${QAS_TYPES}")
  message(DEBUG "QAS_MODEL_DIRNAME: ${QAS_MODEL_DIRNAME}")

  # Set Model library name
  string(TOLOWER "quicc_${QAS_MODEL_DIRNAME}" _model_lib)

  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

  # Set library visibility
  set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

  # Create executables
  foreach(type ${QAS_TYPES})
    set(ModelId "${target}/${type}")
    message(DEBUG "ModelId: ${ModelId}")
    string(TOLOWER "${_model_lib}_${type}" modLib)

    # Config executable
    quicc_setup_exe(${ModelId}
      POSTFIX "StabilityConfig"
      EXESOURCE "${QUICC_TOOLS_DIR}/Stability/WriteConfig.cpp"
      MODELLIB "${modLib}"
      EXTRALIBS QuICC::Tools::Stability)

    # Stability executable
    quicc_setup_exe(${ModelId}
      POSTFIX "Stability"
      EXESOURCE "${QUICC_TOOLS_DIR}/Stability/RunStability.cpp"
      MODELLIB "${modLib}"
      EXTRALIBS QuICC::Tools::Stability)
  endforeach()

  unset(QAS_TYPES)
  unset(modLib)

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()

