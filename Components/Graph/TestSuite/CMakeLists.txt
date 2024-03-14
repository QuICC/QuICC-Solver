include(BundleCatch2)

add_executable(GraphTest
    GraphTests.cpp
)

target_link_libraries(GraphTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
    QuICC::ViewOps
    QuICC::Graph
)

if(TARGET QuICC::Cuda)
    target_link_libraries(GraphTest QuICC::Cuda)
endif()

add_test(
    NAME GraphTest
    COMMAND GraphTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )
