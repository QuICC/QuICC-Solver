include(BundleCatch2)


add_executable(ProfilerTest ProfilerTests.cpp)
target_link_libraries(ProfilerTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Interface
    Catch2::Catch2
    )

add_test(
    NAME ProfilerTest
    COMMAND ProfilerTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )

# Tracker is tested through the profiler interface which may or may not link to the
# Tracker library. In order to always test it, we need to add it manually when missing
if(NOT TARGET ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Tracker)
    message(STATUS "not backend")
    add_subdirectory(../Tracker ${CMAKE_CURRENT_BINARY_DIR}/Tracker EXCLUDE_FROM_ALL)
endif()

add_executable(TrackerTest TrackerTests.cpp)
target_link_libraries(TrackerTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Interface
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Tracker
    Catch2::Catch2
    )

add_test(
    NAME TrackerTest
    COMMAND TrackerTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )


# Instead this test only makes sense if likwid is present
if(TARGET LIKWID::LIKWID)
    add_executable(LikwidTest LikwidTests.cpp)
    target_link_libraries(LikwidTest ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Interface Catch2::Catch2 )

    add_test(
    NAME LikwidTest
    COMMAND likwid-perfctr -C S0:0 -g FLOPS_DP -m LikwidTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"

  )
    if(QUICC_MPI)
        set_tests_properties(LikwidTest PROPERTIES DISABLED True)
    endif()
endif()

# if(QUICC_MPI AND NOT QUICC_MPI_CI)
#     add_executable(ProfilerTestMpi ProfilerTestsMpi.cpp)
#     target_link_libraries(ProfilerTestMpi
#         ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
#         ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}::Tracker
#         MPI::MPI_CXX
#     )

#     # check which command is available
#     foreach(_mpiexe IN ITEMS srun mpirun)
#       message(VERBOSE "_mpiexe: ${_mpiexe}")
#       find_program(mpiexe ${_mpiexe})
#       if(mpiexe STREQUAL "mpiexe-NOTFOUND")
#         message(VERBOSE "not found")
#       else()
#         message(VERBOSE "found")
#         break()
#       endif()
#     endforeach()

#     add_test(
#     NAME ProfilerTestMpi
#     COMMAND ${_mpiexe} -n 4 ProfilerTestMpi
#     WORKING_DIRECTORY
#     "${QUICC_WORK_DIR}"
#   )
# endif()
