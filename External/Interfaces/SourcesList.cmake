# Create list of sources
set(MHDSources
)

if(QUICC_SPLINALG STREQUAL "Pardiso")
   list(APPEND MHDSources Pardiso_Real.cpp Pardiso_Complex.cpp)
endif(QUICC_SPLINALG STREQUAL "Pardiso")
