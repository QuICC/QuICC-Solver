include(BundleCatch2)

add_executable(ProfilerTest ProfilerTests.cpp)
target_link_libraries(ProfilerTest ${QUICC_CURRENT_COMPONENT_LIB}::Interface Catch2::Catch2 )

add_test(
    NAME ProfilerTest
    COMMAND ProfilerTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )

add_executable(TrackerTest TrackerTests.cpp)
target_link_libraries(TrackerTest ${QUICC_CURRENT_COMPONENT_LIB}::Interface Catch2::Catch2 )

# Tracker is tested through the profiler interface which may or may not link to the
# Tracker library. In order to always test it, we need to add it manually when missing
if(NOT TARGET ${QUICC_CURRENT_COMPONENT_LIB}::Tracker)
    add_subdirectory(../Tracker ${CMAKE_CURRENT_BINARY_DIR}/Tracker EXCLUDE_FROM_ALL)
endif()
target_link_libraries(TrackerTest ${QUICC_CURRENT_COMPONENT_LIB}::Tracker)

add_test(
    NAME TrackerTest
    COMMAND TrackerTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )

# Instead this test only makes sense if likwid is present
if(TARGET LIKWID::LIKWID)
    add_executable(LikwidTest LikwidTests.cpp)
    target_link_libraries(LikwidTest ${QUICC_CURRENT_COMPONENT_LIB}::Interface Catch2::Catch2 )

    add_test(
    NAME LikwidTest
    COMMAND likwid-perfctr -C S0:0 -g FLOPS_DP -m LikwidTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )
endif()

# mpi
if(QUICC_MPI)
    add_executable(ProfilerTestMpi ProfilerTestsMpi.cpp)
    target_link_libraries(ProfilerTestMpi ${QUICC_CURRENT_COMPONENT_LIB}::Profiling MPI::MPI_CXX)

    add_test(
    NAME ProfilerTestMpi
    COMMAND mpirun -n 4 ProfilerTestMpi
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )
endif()
