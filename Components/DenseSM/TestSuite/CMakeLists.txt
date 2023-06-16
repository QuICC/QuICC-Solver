include(BundleCatch2)

# Set component test library name
set(QUICC_CURRENT_COMPONENT_TEST_LIB ${QUICC_CURRENT_COMPONENT_LIB}_ts)

# Create component Testsuite library
add_library(${QUICC_CURRENT_COMPONENT_TEST_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_TEST_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(${QUICC_CURRENT_COMPONENT_TEST_LIB}
  ${QUICC_CURRENT_COMPONENT_LIB}
  Catch2::Catch2
)

add_subdirectory(DenseSM)
