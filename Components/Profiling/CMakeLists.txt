# CMake setup for Profiling component
set(QUICC_CURRENT_COMPONENT_LIB Profiling)

add_subdirectory(Timers)
add_subdirectory(Profiler)
if(QUICC_PROFILER_BACKEND STREQUAL "native")
  add_subdirectory(Tracker)
endif()

# Alias Target for all components
add_library(Profiling_Profiling INTERFACE)
target_link_libraries(Profiling_Profiling INTERFACE Profiling::Profiler Profiling::Timers)
add_library(Profiling::Profiling ALIAS Profiling_Profiling)

# Teststuite
option(QUICC_TESTSUITE_PROFILING "Enable Profiling component testsuite?" OFF)
if(QUICC_TESTSUITE_PROFILING)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_PROFILING)
