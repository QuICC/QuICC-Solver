quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ISparseSMOperator.cpp
  IChebyshevOperator.cpp
  IWorlandOperator.cpp
)

add_subdirectory(Chebyshev)
add_subdirectory(Version)
add_subdirectory(Worland)
