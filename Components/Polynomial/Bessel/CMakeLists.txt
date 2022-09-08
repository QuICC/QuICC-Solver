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
  if(NOT BOOST_ROOT)
    message(VERBOSE "setting BOOST_ROOT")
    set(_ALL_PATHS "$ENV{CPATH};$ENV{C_INCLUDE_PATH};$ENV{CPLUS_INCLUDE_PATH}")
    if(_ALL_PATHS STREQUAL "")
      string(REPLACE ":" ";" _ALL_PATHS ${_ALL_PATHS})
      include(ListFindRegex)
      quicc_list(FIND_REGEX _ALL_PATHS "boost" BOOST_ROOT)
      message(VERBOSE "BOOST_ROOT: ${BOOST_ROOT}")
    endif()
  endif()
  find_package(Boost)
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost, try to specify path: -DBOOST_ROOT=</path/to/boost>")
  endif()
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
