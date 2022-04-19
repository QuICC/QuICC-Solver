include(BundleCatch2)

# Set component test library name
set(QUICC_CURRENT_COMPONENT_TEST_LIB ${QUICC_CURRENT_COMPONENT_LIB}_ts)

# Create component Testsuite library
add_library(${QUICC_CURRENT_COMPONENT_TEST_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_TEST_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_TEST_LIB} PUBLIC "include" "Tests")
target_link_libraries(${QUICC_CURRENT_COMPONENT_TEST_LIB}
  ${QUICC_CURRENT_COMPONENT_LIB} Catch2::Catch2
)

# Copy component TestSuite Python files
add_custom_target(${QUICC_CURRENT_COMPONENT_TEST_LIB}_updatepy)
add_custom_command(TARGET ${QUICC_CURRENT_COMPONENT_TEST_LIB}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
   "${CMAKE_CURRENT_SOURCE_DIR}/Python"
   "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite/Python"
  COMMENT "Copying Python files for Common TestSuite"
  VERBATIM
)
add_dependencies(${QUICC_CURRENT_COMPONENT_TEST_LIB} ${QUICC_CURRENT_COMPONENT_TEST_LIB}_updatepy)

# Add source and tests directories
add_subdirectory(src)
add_subdirectory(Tests)
