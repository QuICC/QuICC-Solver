include(BundleCatch2)

add_executable(TypesTest
    TypesTests.cpp
)
target_link_libraries(TypesTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
)

add_test(
    NAME TypesTest
    COMMAND TypesTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
)
