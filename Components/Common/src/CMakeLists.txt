target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Hasher.cpp
  QuICCTimer.cpp
  )

add_subdirectory(Debug)
add_subdirectory(Enums)
add_subdirectory(Timers)
add_subdirectory(Tools)
