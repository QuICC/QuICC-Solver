quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  CouplingFeature.cpp
  CouplingInformation.cpp
  CorrectSolution.cpp
  Dispatchers.cpp
  IEquation.cpp
  IFieldEquation.cpp
  IScalarEquation.cpp
  IVectorEquation.cpp
  EquationData.cpp
  EquationOptions.cpp
  EquationParameters.cpp
  )

add_subdirectory(details)
add_subdirectory(Tools)
