# CMake setup for Framework component
set(QUICC_CURRENT_COMPONENT_DIR Components/Framework)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_FRAMEWORK_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

include(cmake.d/setup/RegisterTransformForward.cmake)
include(cmake.d/setup/RegisterTransformBackward.cmake)
include(cmake.d/setup/RegisterTransformReductor.cmake)
include(cmake.d/setup/RegisterTransformPath.cmake)
include(cmake.d/setup/RegisterTimestepId.cmake)
include(cmake.d/setup/RegisterTimestepRegister.cmake)
include(cmake.d/setup/RegisterBcScheme.cmake)
include(cmake.d/setup/RegisterBcName.cmake)

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
  "SHELL:${CMAKE_CXX_FLAGS}")

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
  ${QUICC_POLYNOMIAL_COMPONENT_LIB}
  ${QUICC_PYQUICC_COMPONENT_LIB}
  QuICC::Io
  QuICC::Common
  QuICC::Framework::GitHash
  QuICC::Profiling
  )

# Transform component
add_library(${QUICC_CURRENT_COMPONENT_LIB}_transform "")
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_transform PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_transform PUBLIC QuICC::Transform)

# Link to main lib
target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC ${QUICC_CURRENT_COMPONENT_LIB}_transform)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}_transform
  COMPONENT FrameworkTransform
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    ${PROJECT_BINARY_DIR}/include/QuICC
    )

# Add source directory
add_subdirectory(src)

# Kokkos
if(QUICC_USE_KOKKOS)
  include(setup/KokkosAddLibrary)
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_transform PRIVATE QuICC::Kokkos)
  # fix file properties
  quicc_fix_target_sources_kokkos(${QUICC_CURRENT_COMPONENT_LIB}_transform)
endif()

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Alias
add_library(${QUICC_NAMESPACE}Framework
  ALIAS ${QUICC_CURRENT_COMPONENT_LIB})

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Framework
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    ${PROJECT_BINARY_DIR}/include/QuICC
)

# Export targets
install(EXPORT QuICCTargets
  FILE QuICCTargets.cmake
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
