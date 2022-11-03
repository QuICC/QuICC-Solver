quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Hasher.cpp
  Precision.cpp
  QuICCEnv.cpp
  QuICCTimer.cpp
  )

add_subdirectory(Debug)
add_subdirectory(Enums)
add_subdirectory(Environment)
add_subdirectory(Framework)
add_subdirectory(Timers)
add_subdirectory(Tools)
add_subdirectory(Version)
