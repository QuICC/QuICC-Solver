###################################################
#------------------ Components -------------------#
###################################################

message(STATUS "Components setup")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_NAMESPACE "QuICC::")

if(NOT QUICC_USE_SYSTEM_QUICC)
    add_subdirectory("Components")
else()
    message(STATUS "nothing to setup, using installed framework")
endif()


list(POP_BACK CMAKE_MESSAGE_INDENT)
