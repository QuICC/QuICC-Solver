# CMake setup for Profiling component
set(QUICC_CURRENT_COMPONENT_LIB Profiling)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

add_subdirectory(Timers)
add_subdirectory(Profiler)
if(QUICC_PROFILE_BACKEND STREQUAL "native")
  add_subdirectory(Tracker)
endif()

# Inteface target for all components
add_library("${QUICC_CURRENT_COMPONENT_LIB}" INTERFACE)
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}"
  INTERFACE
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}Interface
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}Timers)
# Backends
if(QUICC_PROFILE_BACKEND STREQUAL "native")
  message(VERBOSE "profiler using native backend")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}"
    INTERFACE ${QUICC_CURRENT_COMPONENT_LIB}::Tracker)
elseif(QUICC_PROFILE_BACKEND STREQUAL "likwid")
  message(VERBOSE "profiler using likwid backend")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}"
    INTERFACE LIKWID::LIKWID)
endif()
# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
option(QUICC_TESTSUITE_PROFILING "Enable Profiling component testsuite?" OFF)
if(QUICC_TESTSUITE_PROFILING)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_PROFILING)


# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
