# CMake setup for TimestepPc
set(QUICC_CURRENT_COMPONENT_LIB FiniteDiff)
set(QUICC_CURRENT_COMPONENT_DIR Components/FiniteDiff)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Components/Polynomial/cmake.d")

add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

quicc_create_option(NAME QUICC_FINITEDIFF_SPHERE_ORDER
                    OPTS "2" "4" "6" "8"
                    LABEL "Order of sphere finite differences"
                    ADVANCED)

# Default finite different order
target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC QUICC_FINITEDIFF_SPHERE_ORDER=${QUICC_FINITEDIFF_SPHERE_ORDER}
  )

set(PHS
  FiniteDiff/Sphere/UniformRadialGrid.hpp
  FiniteDiff/Sphere/Operator.hpp
  FiniteDiff/Sphere/Id.hpp
  FiniteDiff/Sphere/R1.hpp
  FiniteDiff/Sphere/D1.hpp
  FiniteDiff/Sphere/Overr1.hpp
  FiniteDiff/Sphere/Overr1D1R1.hpp
  FiniteDiff/Sphere/SLapl.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

add_subdirectory(FiniteDiff)

# Dependencies
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    QuICC::Types
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# TestSuite
include(setup/TestSuite)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
)
