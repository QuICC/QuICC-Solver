quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Setup.cpp
  Tools.cpp
)

add_subdirectory(ALegendre)
add_subdirectory(Worland)
