#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

message(STATUS "Framework setup")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")


###################################################
#----------------- MEMORY USAGE ------------------#
###################################################

#
# Choose the type of memory usage setup the code is using.
# Possible options are: High, Low
#
quicc_create_option(NAME QUICC_MEMORYUSAGE
                    OPTS "High" "Low"
                    LABEL "Memory usage")
quicc_add_definition(QUICC_MEMORYUSAGE)


###################################################
#------------------- PROFILING -------------------#
###################################################
quicc_create_option(NAME QUICC_PROFILE_BACKEND
  OPTS none native likwid
  LABEL "Profiler backend."
  ADVANCED)

if(QUICC_PROFILE_BACKEND STREQUAL "likwid")
  find_package(LIKWID REQUIRED)
endif()

if(QUICC_PROFILE_BACKEND STREQUAL "native")
  quicc_create_option(NAME QUICC_PROFILE_NATIVE_WRITER
    OPTS HighFive none
    LABEL "Native profiler writer backend."
    ADVANCED)
  if(QUICC_PROFILE_NATIVE_WRITER STREQUAL "HighFive")
      include(BundleHighFive)
  endif()
  set(QUICC_PROFILE_NATIVE_SAMPLE 100 CACHE STRING "Native profiler sample size.")
  mark_as_advanced(QUICC_PROFILE_NATIVE_WRITER)
else()
  unset(QUICC_PROFILE_NATIVE_WRITER CACHE)
endif()

quicc_create_option(NAME QUICC_PROFILE_LEVEL
  OPTS 0 1 2 3 4
  LABEL "Profiler granularity."
  ADVANCED)

list(POP_BACK CMAKE_MESSAGE_INDENT)

