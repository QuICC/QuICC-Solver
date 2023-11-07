target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IWorlandProjector.cpp
)

add_subdirectory(Base)
add_subdirectory(Kokkos)
add_subdirectory(KokkosCuda)
