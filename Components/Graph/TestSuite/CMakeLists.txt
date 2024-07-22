include(BundleCatch2)

add_executable(GraphTests
    GraphTestsMain.cpp
    GraphTests.cpp
)

target_link_libraries(GraphTests
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    QuICC::Graph
    QuICC::TestSuite
    Catch2::Catch2
)

if(TARGET QuICC::Cuda)
    target_sources(GraphTests
        PRIVATE
            GraphTestsCuda.cpp
    )
    target_link_libraries(GraphTests QuICC::Cuda)
endif()

add_test(
    NAME GraphTests
    COMMAND GraphTests
    WORKING_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )

# if(QUICC_USE_MPI)
    add_executable(GraphTestsMpi
        GraphTestsMpi.cpp
    )

    target_link_libraries(GraphTestsMpi
        ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
        QuICC::Graph
        QuICC::Environment
        QuICC::TestSuite
        Catch2::Catch2
    )

    add_test(
        NAME GraphTestsMpi
        COMMAND GraphTestsMpi
        WORKING_DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
# endif()