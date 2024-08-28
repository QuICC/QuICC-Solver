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

if(QUICC_USE_MPI)
    add_executable(GraphTestsMpi
        GraphTestsMpi.cpp
    )

    target_link_libraries(GraphTestsMpi
        ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
        QuICC::Graph
        QuICC::Environment
        QuICC::TestSuite
        QuICC::Profiling
        Catch2::Catch2
    )

    # Reusing LoadSplitter data
    set(QUICC_WORK_DIR "${CMAKE_BINARY_DIR}/Components/Framework/TestSuite")

    # Base dir for reference data archive
    set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

    # Fetch reference data
    include(FetchTestReference)
    quicc_fetch_test_reference(
        ViewTransposeMpiTest
        COMPONENT Framework
        FILENAME "LoadSplitter.tar.gz"
        ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
        DATADIR ${QUICC_WORK_DIR}
    )

    # Copy graph files
    file(GLOB mlirInputs
        "${CMAKE_CURRENT_SOURCE_DIR}/*.mlir"
    )
    file(COPY ${mlirInputs}
        DESTINATION "${QUICC_WORK_DIR}")

    add_test(
        NAME GraphTestsMpi
        COMMAND GraphTestsMpi
        WORKING_DIRECTORY
        "${QUICC_WORK_DIR}"
    )
endif()
