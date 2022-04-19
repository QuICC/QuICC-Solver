###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

#
# General compiler Settings
#

# TODO, add compiler detection and add specific flags here, for instance -march=native, -ggdb3

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -DEIGEN_NO_DEBUG -DNDEBUG -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -DQUICC_DEBUG" CACHE STRING "" FORCE)
if(QUICC_DISABLE_RDYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(QUICC_DISABLE_RDYNAMIC)
if(QUICC_ENABLE_DYNAMIC)
   set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-dynamic ${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}" CACHE STRING "" FORCE)
   mark_as_advanced(FORCE CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(QUICC_ENABLE_DYNAMIC)

# Threads
if(${QUICC_THREADS} STREQUAL "pthread")
   add_link_options(-pthread)
elseif(${QUICC_THREADS} STREQUAL "openmp")
   add_link_options(-fopenmp)
endif()
