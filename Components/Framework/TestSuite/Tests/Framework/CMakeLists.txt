list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
message(VERBOSE "TestSuite")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
add_subdirectory(Communicators)
add_subdirectory(LoadSplitter)
add_subdirectory(TransformConfigurators)
add_subdirectory(TransformCoordinator)
add_subdirectory(Transform)
list(POP_BACK CMAKE_MESSAGE_INDENT)

list(POP_BACK CMAKE_MESSAGE_INDENT)
