# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB}_SphericalBessel "Jnl.cpp")
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_SphericalBessel
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

add_library(Polynomial::SphericalBessel ALIAS
  ${QUICC_CURRENT_COMPONENT_LIB}_SphericalBessel)

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_SphericalBessel PUBLIC
  QuICC::Types
)

# Teststuite
#if(QUICC_TESTSUITE_POLYNOMIAL)
#  add_subdirectory(TestSuite)
#endif(QUICC_TESTSUITE_POLYNOMIAL)

# Export
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}_SphericalBessel
  COMPONENT SphericalBessel
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
)
