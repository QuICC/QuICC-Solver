quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Hasher.cpp
  Precision.cpp
  QuICCEnv.cpp
  QuICCTimer.cpp
  )

add_subdirectory(Arithmetics)
add_subdirectory(Debug)
add_subdirectory(Enums)
add_subdirectory(Environment)
add_subdirectory(Framework)
add_subdirectory(ModelOperator)
add_subdirectory(ModelOperatorBoundary)
add_subdirectory(NonDimensional)
add_subdirectory(PhysicalNames)
add_subdirectory(PseudospectralTag)
add_subdirectory(RuntimeStatus)
add_subdirectory(SolveTiming)
add_subdirectory(Timers)
add_subdirectory(Tools)
add_subdirectory(Version)
