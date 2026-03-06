include(BundleCatch2)

# Set component test library name
set(QUICC_CURRENT_COMPONENT_TEST_LIB ${QUICC_CURRENT_COMPONENT_LIB}_ts)

# Create component Testsuite library
add_library(${QUICC_CURRENT_COMPONENT_TEST_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_TEST_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_TEST_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(${QUICC_CURRENT_COMPONENT_TEST_LIB}
  PUBLIC
    ${QUICC_CURRENT_COMPONENT_LIB}
    QuICC::TestSuite
    Catch2::Catch2
)

# Add source and tests directories
add_subdirectory(Timestep)
add_subdirectory(Tests)
