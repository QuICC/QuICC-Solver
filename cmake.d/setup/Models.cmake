###################################################
#-------------------- Models ---------------------#
###################################################

message(STATUS "Models setup")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

if(EXISTS "${QUICC_MODELS_DIR}/CMakeLists.txt")
  add_subdirectory(Models)
endif(EXISTS "${QUICC_MODELS_DIR}/CMakeLists.txt")

list(POP_BACK CMAKE_MESSAGE_INDENT)
