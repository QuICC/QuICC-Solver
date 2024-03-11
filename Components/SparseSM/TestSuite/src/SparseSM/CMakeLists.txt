quicc_target_sources(${QUICC_CURRENT_COMPONENT_TEST_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Io.cpp
  TestArgs.cpp
  TestUtils.cpp
)

add_subdirectory(Bessel)
add_subdirectory(Chebyshev)
add_subdirectory(Worland)
