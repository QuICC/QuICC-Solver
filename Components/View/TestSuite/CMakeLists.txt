include(BundleCatch2)


add_executable(ViewTest 
    ViewTests.cpp
    ViewDenseTests.cpp 
    ViewSparseTests.cpp 
    ViewTests.cu
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
