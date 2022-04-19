include(BundleCatch2)

find_package(Boost REQUIRED)

add_executable(BesselTest BesselTests.cpp)
target_link_libraries(BesselTest Polynomial::Bessel Catch2::Catch2 Boost::headers)

add_test(
    NAME BesselTest
    COMMAND BesselTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )

