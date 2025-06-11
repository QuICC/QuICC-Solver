message(VERBOSE "TestSuite")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

add_subdirectory(Chebyshev)
add_subdirectory(ALegendre)
add_subdirectory(Fourier)
add_subdirectory(Worland)

list(POP_BACK CMAKE_MESSAGE_INDENT)
