quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ThreeTermRecurrence.cpp
)

add_subdirectory(ALegendre)
add_subdirectory(Jacobi)
add_subdirectory(Legendre)
add_subdirectory(Quadrature)
add_subdirectory(Version)
add_subdirectory(Worland)
