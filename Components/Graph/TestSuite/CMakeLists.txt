include(BundleCatch2)

add_executable(GraphTest
    GraphTestsMain.cpp
    GraphTests.cpp
)

target_link_libraries(GraphTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
    QuICC::Graph
)

if(TARGET QuICC::Cuda)
    target_sources(GraphTest
        PRIVATE
            GraphTestsCuda.cpp
    )
    target_link_libraries(GraphTest QuICC::Cuda)
endif()

add_test(
    NAME GraphTest
    COMMAND GraphTest
    WORKING_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )
