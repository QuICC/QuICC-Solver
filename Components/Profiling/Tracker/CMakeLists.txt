# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB}_Tracker STATIC
  Tracker.cpp
  MemoryTracker.cpp)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_Tracker PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_Tracker PRIVATE ${QUICC_CURRENT_COMPONENT_LIB}::Timers)

# Add MPI dependency
if(QUICC_MPI)
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_Tracker PRIVATE MPI::MPI_CXX)
endif(QUICC_MPI)

add_library(${QUICC_CURRENT_COMPONENT_LIB}::Tracker ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_Tracker)
