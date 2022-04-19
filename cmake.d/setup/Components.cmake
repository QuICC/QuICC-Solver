###################################################
#------------------ Components -------------------#
###################################################

message(STATUS "Components setup")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

add_subdirectory("Components")

list(POP_BACK CMAKE_MESSAGE_INDENT)
