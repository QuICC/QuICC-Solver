# CMake setup for Types
set(QUICC_CURRENT_COMPONENT_LIB Types)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

add_library(${QUICC_CURRENT_COMPONENT_LIB}
  STATIC
    Types/Precision.cpp
  )

# Options
include(cmake.d/Options.cmake)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set(PHS
  Types/BasicTypes.hpp
  Types/Constants.hpp
  Types/DecoupledComplex.hpp
  Types/DecoupledComplexInternal.hpp
  Types/MatrixOperationsInternal.hpp
  Types/MpTypedefs.hpp
  Types/Precision.hpp
  Types/Typedefs.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Dependencies
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
      Eigen3::Eigen
)

# Multiple precision def
# target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}
#     INTERFACE QUICC_HAS_CUDA_BACKEND)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# # Teststuite
# option(QUICC_TESTSUITE_Types "Enable Types component testsuite?" OFF)
# if(QUICC_TESTSUITE_Types)
#   add_subdirectory(TestSuite)
# endif(QUICC_TESTSUITE_Types)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
)
