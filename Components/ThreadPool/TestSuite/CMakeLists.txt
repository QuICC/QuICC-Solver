include(BundleCatch2)

add_executable(ThreadPoolTests
    ThreadPoolTests.cpp
    ThreadPoolTestsMain.cpp
)

target_link_libraries(ThreadPoolTests
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    QuICC::TestSuite
    QuICC::Profiling
    Catch2::Catch2
)

add_test(
    NAME ThreadPoolTests
    COMMAND ThreadPoolTests
    WORKING_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )
