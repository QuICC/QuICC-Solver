include(BundleCatch2)


add_executable(ViewTest
    ViewTests.cpp
    ViewTests.cu
    ViewDenseTests.cpp
    ViewSparseTests.cpp
    ViewMapsTests.cpp
)
target_link_libraries(ViewTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
    QuICC::Memory
    )

if(TARGET QuICC::Cuda)
    target_link_libraries(ViewTest QuICC::Cuda)
endif()

add_test(
    NAME ViewTest
    COMMAND ViewTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
  )
