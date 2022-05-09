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
   add_compile_definitions("QUICC_PROFILE")

   if(QUICC_MPI)
      option(QUICC_PROFILE_PERCORE "Write per core profiling data?" OFF)
      if(QUICC_PROFILE_PERCORE)
         add_compile_definitions("QUICC_PROFILE_PERCORE")
      endif(QUICC_PROFILE_PERCORE)
   endif(QUICC_MPI)

endif(QUICC_PROFILE)

###################################################
#--------------- STORAGE PROFILING ---------------#
###################################################

#
# Used storage requirements profiler?
#
option(QUICC_STORAGEPROFILE "Activate internal storage profiler?" OFF)

if(QUICC_STORAGEPROFILE)
   add_compile_definitions("QUICC_STORAGEPROFILE")
endif(QUICC_STORAGEPROFILE)

list(POP_BACK CMAKE_MESSAGE_INDENT)
