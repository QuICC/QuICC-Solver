target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  Utils.cpp
  ThreadPoolBoost.cpp
  QuICCThreads.cpp
  )

add_subdirectory(details)
