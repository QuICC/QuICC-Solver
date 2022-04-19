message(VERBOSE "TestSuite")

list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
add_subdirectory(Quadrature)
add_subdirectory(ALegendre)
add_subdirectory(Jacobi)
add_subdirectory(Legendre)
add_subdirectory(Worland)
list(POP_BACK CMAKE_MESSAGE_INDENT)
