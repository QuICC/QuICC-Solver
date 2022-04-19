# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB}_Interface INTERFACE)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_Interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Add MPI dependency
if(QUICC_MPI)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Interface" INTERFACE MPI::MPI_CXX)
endif(QUICC_MPI)

# Backends
if(QUICC_PROFILER_BACKEND STREQUAL "native")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Interface" INTERFACE ${QUICC_CURRENT_COMPONENT_LIB}::Tracker)
elseif(QUICC_PROFILER_BACKEND STREQUAL "likwid")
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}_Interface" INTERFACE LIKWID::LIKWID)
endif()

add_library(${QUICC_CURRENT_COMPONENT_LIB}::Interface ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_Interface)
