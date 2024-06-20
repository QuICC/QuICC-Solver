# CMake setup for Environment (serial or parallel)
set(QUICC_CURRENT_COMPONENT_LIB Environment)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )


add_library(${QUICC_CURRENT_COMPONENT_LIB} STATIC "")

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

add_subdirectory(Environment)

# Add MPI dependency
if(QUICC_USE_MPI)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}"
    PUBLIC
      MPI::MPI_CXX
  )
  # to do, change to QUICC_USE_MPI
  target_compile_definitions("${QUICC_CURRENT_COMPONENT_LIB}"
    PUBLIC
      QUICC_MPI
  )
endif()

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
  ALIAS
    "${QUICC_CURRENT_COMPONENT_LIB}")

# # Teststuite
# option(QUICC_TESTSUITE_ENVIRONMENT "Enable Environment component testsuite?" OFF)
# if(QUICC_TESTSUITE_ENVIRONMENT)
#   add_subdirectory(TestSuite)
# endif(QUICC_TESTSUITE_ENVIRONMENT)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
