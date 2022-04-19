quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  CouplingFeature.cpp
  CouplingInformation.cpp
  IEquation.cpp
  IFieldEquation.cpp
  IScalarEquation.cpp
  IVectorEquation.cpp
  EquationData.cpp
  EquationParameters.cpp
  )

add_subdirectory(Tools)
