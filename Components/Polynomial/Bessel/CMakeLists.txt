quicc_create_option(NAME QUICC_BESSEL_IMPL
                    OPTS "Boost" "std"
                    LABEL "Bessel functions backend"
                    ADVANCED)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB}_Bessel INTERFACE)
# target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_Bessel INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )
quicc_target_add_definition(${QUICC_CURRENT_COMPONENT_LIB}_Bessel INTERFACE
  OPTION QUICC_BESSEL_IMPL)

if(QUICC_BESSEL_IMPL STREQUAL "Boost")
  include(QuICCFindBoostRequired)
  quicc_find_boost_required()
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_Bessel INTERFACE Boost::headers)
endif()

add_library(Polynomial::Bessel ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_Bessel)

# Teststuite
if(QUICC_TESTSUITE_POLYNOMIAL)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_POLYNOMIAL)

# Export
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}_Bessel
  COMPONENT Bessel
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
)
