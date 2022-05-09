# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB}_Interface INTERFACE)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_Interface
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Add MPI dependency
if(QUICC_MPI)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Interface"
    INTERFACE MPI::MPI_CXX)
endif(QUICC_MPI)

# Granularity level
target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}_Interface
  INTERFACE QUICC_PROFILE_LEVEL=${QUICC_PROFILE_LEVEL}
  )

# Backend
quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB}_Interface INTERFACE
  OPTION QUICC_PROFILE_BACKEND)

# Alias
add_library(${QUICC_CURRENT_COMPONENT_LIB}::Interface
  ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_Interface)

