#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

###################################################
#-------------------- DEBUGGER -------------------#
###################################################

message(STATUS "Debug setup")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
   set(QUICC_DEBUG ON)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if(QUICC_DEBUG)
   message(VERBOSE "Debug: Active")
endif(QUICC_DEBUG)

if(QUICC_DEBUG)
   add_definitions("-DQUICC_DEBUG")
else(QUICC_DEBUG)
   add_definitions("-DQUICC_NO_DEBUG")
endif(QUICC_DEBUG)

###################################################
#------------------- PROFILING -------------------#
###################################################

#
# Use internal profiler and type of profiler
#
option(QUICC_PROFILE "Activate internal profiler?" OFF)


if(QUICC_PROFILE)
   add_definitions("-DQUICC_PROFILE")

   if(QUICC_MPI)
      option(QUICC_PROFILE_PERCORE "Write per core profiling data?" OFF)
      if(QUICC_PROFILE_PERCORE)
         add_definitions("-DQUICC_PROFILE_PERCORE")
      endif(QUICC_PROFILE_PERCORE)
   endif(QUICC_MPI)

   quicc_create_option(NAME QUICC_PROFILER_BACKEND
                    OPTS native likwid
                    LABEL "Profiler backend."
                    ADVANCED)
   quicc_add_definition(QUICC_PROFILER_BACKEND)

   if(QUICC_PROFILER_BACKEND STREQUAL "likwid")
      find_package(LIKWID REQUIRED)
   endif()

   quicc_create_option(NAME QUICC_PROFILER_LEVEL
                    OPTS 0 1 2 3
                    LABEL "Profiler granularity."
                    ADVANCED)
   quicc_add_definition(QUICC_PROFILER_LEVEL)


endif(QUICC_PROFILE)

###################################################
#--------------- STORAGE PROFILING ---------------#
###################################################

#
# Used storage requirements profiler?
#
option(QUICC_STORAGEPROFILE "Activate internal storage profiler?" OFF)

if(QUICC_STORAGEPROFILE)
   add_definitions("-DQUICC_STORAGEPROFILE")
endif(QUICC_STORAGEPROFILE)

list(POP_BACK CMAKE_MESSAGE_INDENT)
