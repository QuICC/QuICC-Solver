# CMake setup for Framework component
set(QUICC_CURRENT_COMPONENT_DIR Components/Framework)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_FRAMEWORK_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

include(cmake.d/setup/RegisterTransformForward.cmake)
include(cmake.d/setup/RegisterTransformBackward.cmake)
include(cmake.d/setup/RegisterTransformReductor.cmake)
include(cmake.d/setup/RegisterTransformPath.cmake)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Add fancy flags as requirement for the library
target_compile_options(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${CMAKE_CXX_FLAGS})

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_TRANSFORM_COMPONENT_LIB}
  ${QUICC_POLYNOMIAL_COMPONENT_LIB}
  ${QUICC_PYQUICC_COMPONENT_LIB}
  QuICC::Io
  QuICC::Common
  QuICC::Profiling
  )

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Framework
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    ${PROJECT_BINARY_DIR}/include/QuICC
)

# Export targets
install(EXPORT QuICCExpTargets
  FILE QuICCExpTargets.cmake
  NAMESPACE ${QUICC_NAMESPACE}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/quicc/cmake
)

# Create config file
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/quiccConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/share/quicc/cmake
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/quiccConfigVersion.cmake
  VERSION 0.0.0
  COMPATIBILITY SameMajorVersion )

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/quiccConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/quiccConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/quicc/cmake
)

# Add library alias (as the imported one)
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
  ALIAS ${QUICC_CURRENT_COMPONENT_LIB})
