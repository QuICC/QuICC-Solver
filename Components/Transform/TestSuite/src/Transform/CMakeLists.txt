quicc_target_sources(${QUICC_CURRENT_COMPONENT_TEST_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  TestArgs.cpp
  Tools.cpp
)

add_subdirectory(Chebyshev)
add_subdirectory(ALegendre)
add_subdirectory(Fourier)
add_subdirectory(Worland)
