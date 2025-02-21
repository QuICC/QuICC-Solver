###################################################
#--------------------- Tools ---------------------#
###################################################

message(STATUS "Tools setup")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

if(EXISTS "${QUICC_TOOLS_DIR}/CMakeLists.txt")
  add_subdirectory(Tools)
endif(EXISTS "${QUICC_TOOLS_DIR}/CMakeLists.txt")

list(POP_BACK CMAKE_MESSAGE_INDENT)
