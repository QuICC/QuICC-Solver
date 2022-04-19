add_library(${QUICC_CURRENT_COMPONENT_LIB}_Timers STATIC
  TimerTools.cpp
  SteadyTimer.cpp
  )

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_Timers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Add in the parallel or serial timers depending on setup
if(QUICC_MPI)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_Timers PRIVATE
    MpiTimer.cpp
    )
else(QUICC_MPI)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_Timers PRIVATE
    SerialTimer.cpp
    )
endif(QUICC_MPI)


add_library(${QUICC_CURRENT_COMPONENT_LIB}::Timers ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_Timers)
