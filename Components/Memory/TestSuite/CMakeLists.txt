include(BundleCatch2)

add_executable(MemoryTest MemoryTests.cpp MemoryTests.cu)
target_link_libraries(MemoryTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
)

if(TARGET QuICC::Cuda)
    target_link_libraries(MemoryTest QuICC::Cuda)
endif()

add_test(
    NAME MemoryTest
    COMMAND MemoryTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
)
