quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IDiags.cpp
  Id.cpp
  IdDiags.cpp
  SphLapl.cpp
  SphLaplDiags.cpp
  SphLapl2.cpp
  SphLapl2Diags.cpp
)

add_subdirectory(Value)
add_subdirectory(Insulating)
