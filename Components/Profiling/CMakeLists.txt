# Setup
quicc_create_option(NAME QUICC_PROFILE_BACKEND
                OPTS none native likwid
                LABEL "Profiler backend."
                ADVANCED)

if(QUICC_PROFILE_BACKEND STREQUAL "likwid")
  find_package(LIKWID REQUIRED)
endif()

quicc_create_option(NAME QUICC_PROFILE_LEVEL
                OPTS 0 1 2 3 4
                LABEL "Profiler granularity."
                ADVANCED)


# CMake setup for Profiling component
set(QUICC_CURRENT_COMPONENT_LIB Profiling)

add_subdirectory(Timers)
add_subdirectory(Profiler)
if(QUICC_PROFILE_BACKEND STREQUAL "native")
  add_subdirectory(Tracker)
endif()

# Inteface target for all components
add_library("${QUICC_CURRENT_COMPONENT_LIB}_Profiling" INTERFACE)
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Profiling" INTERFACE Profiling::Interface Profiling::Timers)
# Backends
if(QUICC_PROFILE_BACKEND STREQUAL "native")
  message(VERBOSE "profiler using native backend")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Profiling" INTERFACE ${QUICC_CURRENT_COMPONENT_LIB}::Tracker)
elseif(QUICC_PROFILE_BACKEND STREQUAL "likwid")
  message(VERBOSE "profiler using likwid backend")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Profiling" INTERFACE LIKWID::LIKWID)
endif()
# Alias
add_library(Profiling::Profiling ALIAS "${QUICC_CURRENT_COMPONENT_LIB}_Profiling")

# Teststuite
option(QUICC_TESTSUITE_PROFILING "Enable Profiling component testsuite?" OFF)
if(QUICC_TESTSUITE_PROFILING)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_PROFILING)
