###################################################
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Enable MPI parallelisation
#
option(QUICC_USE_MPI "Use MPI algorithm" OFF)

if(QUICC_USE_MPI)
   if(NOT MPI_FOUND)
      message(SEND_ERROR "QUICC_USE_MPI requires MPI.")
   endif(NOT MPI_FOUND)
   target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC QUICC_MPI)
   option(QUICC_MPI_CI "Running benchmark on the CI?" OFF)
   mark_as_advanced(QUICC_MPI_CI)
endif()


###################################################
#--------------- MPI DATA PACKING ----------------#
###################################################

if(QUICC_USE_MPI)
    quicc_create_option(NAME QUICC_MPIPACK
                        OPTS "MPI" "Manual"
                        LABEL "MPI data packing"
                        ADVANCED)
    quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
    OPTION QUICC_MPIPACK)
endif()


###################################################
#--------------- MPI COMMUNICATION ---------------#
###################################################

if(QUICC_USE_MPI)
    quicc_create_option(NAME QUICC_MPICOMM
                        OPTS "AllToAll" "SendRecv"
                        LABEL "MPI communication"
                        ADVANCED)
    quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
    OPTION QUICC_MPICOMM)
endif()

