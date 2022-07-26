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
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Choose the type of MPI parallelisation or serial setup.
# Possible options are: Serial, Single1D, Single2D, Tubular
#
quicc_create_option(NAME QUICC_MPIALGO
                    OPTS "Serial" "Tubular" "Single1D" "Single2D" "Coupled2D" "SerialMpi"
                    LABEL "MPI algorithm")

if(QUICC_MPIALGO STREQUAL "SerialMpi")
   # Overwrite def
   quicc_add_definition(QUICC_MPIALGO FORCE "Serial")
else()
   quicc_add_definition(QUICC_MPIALGO)
endif()

if(NOT QUICC_MPIALGO STREQUAL "Serial")
   if(NOT MPI_FOUND)
      message(SEND_ERROR "QUICC_MPIALGO ${QUICC_MPIALGO} requires MPI.")
   endif(NOT MPI_FOUND)
   set(QUICC_MPI ON)
   add_definitions("-DQUICC_MPI")
   option(QUICC_MPI_CI "Running benchmark on the CI?" OFF)
   mark_as_advanced(QUICC_MPI_CI)
else()
   set(QUICC_MPI OFF)
endif(NOT QUICC_MPIALGO STREQUAL "Serial")

###################################################
#---------- MPI COMMUNICATION GROUPING -----------#
###################################################

#
# Choose the type of Serial/MPI parallelisation grouping setup.
# Possible options are: Serial, Single1D, Single2D, Tubular
#
set(QUICC_GROUPERS_SERIAL "Equation")
set(QUICC_GROUPERS_SINGLE1D "Single1D" "Equation")
set(QUICC_GROUPERS_SINGLE2D "Single2D" "Equation")
set(QUICC_GROUPERS_TUBULAR "Transform" "Equation" "Single1D" "Single2D")
set(QUICC_GROUPERS_COUPLED2D "Equation" "Single1D")
string(TOUPPER "QUICC_GROUPERS_${QUICC_MPIALGO}" upGrouper)

quicc_create_option(NAME QUICC_TRANSGROUPER
                    OPTS ${${upGrouper}}
                    LABEL "MPI algorithm")
quicc_add_definition(QUICC_TRANSGROUPER)

###################################################
#--------------- MPI DATA PACKING ----------------#
###################################################

if(QUICC_MPI)
   quicc_create_option(NAME QUICC_MPIPACK
                       OPTS "MPI" "Manual"
                       LABEL "MPI data packing"
                       ADVANCED)
   quicc_add_definition(QUICC_MPIPACK)
endif(QUICC_MPI)


###################################################
#--------------- MPI COMMUNICATION ---------------#
###################################################

if(QUICC_MPI)
   quicc_create_option(NAME QUICC_MPICOMM
                       OPTS "AllToAll" "SendRecv"
                       LABEL "MPI communication"
                       ADVANCED)
   quicc_add_definition(QUICC_MPICOMM)
endif(QUICC_MPI)


###################################################
#--------------- TIME INTEGRATORS ----------------#
###################################################

#
# Choose the type of time integration the code is using.
# Possible options are: ImExRKCB2, ImExRKCB3a, ImExRKCB3b, ImExRKCB3c, ImExRKCB3d, ImExRKCB3e, ImExRKCB3f, ImExRKCB4, ImExRK3, ImExSBDF2
#
quicc_create_option(NAME QUICC_TIMESTEPPER
                    OPTS "ImExRKCB2" "ImExRKCB3a" "ImExRKCB3b" "ImExRKCB3c" "ImExRKCB3d" "ImExRKCB3e" "ImExRKCB3f" "ImExRKCB4" "ImExRK3" "ImExSBDF2"
                    LABEL "Time integrator")
quicc_add_definition(QUICC_TIMESTEPPER)




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
endif()

quicc_create_option(NAME QUICC_PROFILE_LEVEL
  OPTS 0 1 2 3 4
  LABEL "Profiler granularity."
  ADVANCED)

list(POP_BACK CMAKE_MESSAGE_INDENT)

