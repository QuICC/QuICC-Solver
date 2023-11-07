target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IWorlandReductor.cpp
  IWorlandPower.cpp
  IWorlandRadialPower.cpp
)

add_subdirectory(Base)
add_subdirectory(Kokkos)
add_subdirectory(KokkosCuda)
