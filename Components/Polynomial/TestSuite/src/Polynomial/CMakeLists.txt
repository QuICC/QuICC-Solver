quicc_target_sources(${QUICC_CURRENT_COMPONENT_TEST_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Io.cpp
  TestArgs.cpp
)

add_subdirectory(Quadrature)
add_subdirectory(ALegendre)
add_subdirectory(Jacobi)
add_subdirectory(Legendre)
add_subdirectory(Worland)
add_subdirectory(Bessel)
