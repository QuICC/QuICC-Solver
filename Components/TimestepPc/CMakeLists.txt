# CMake setup for TimestepPc
set(QUICC_CURRENT_COMPONENT_LIB TimestepPc)
set(QUICC_CURRENT_COMPONENT_DIR Components/TimestepPc)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

include(cmake.d/setup/RegisterTimestepId.cmake)

add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set(PHS
  Timestep/PredictorCorrector/Factory.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

add_subdirectory(Timestep)

# Dependencies
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    QuICC::Framework
)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# TestSuite
string(TOUPPER "${QUICC_CURRENT_COMPONENT_LIB}" _tstag) 
option(QUICC_TESTSUITE_${_tstag} "Enable Timestep/PredictorCorrector component testsuite?" OFF)
if(QUICC_TESTSUITE_${_tstag})
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_${_tstag})

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
)
