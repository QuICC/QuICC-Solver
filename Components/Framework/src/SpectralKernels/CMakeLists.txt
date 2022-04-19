quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  ISpectralKernel.cpp
  IoValidation.cpp
  MakeConstant.cpp
  MakeRandom.cpp
  Set3DModes.cpp
)

add_subdirectory(Sphere)
