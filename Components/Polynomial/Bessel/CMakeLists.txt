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

# Export info
install(TARGETS ${QUICC_CURRENT_COMPONENT_LIB}_Bessel
        EXPORT ${QUICC_CURRENT_COMPONENT_LIB}_BesselTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install header files
# install(FILES MathFunctions.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ${QUICC_CURRENT_COMPONENT_LIB}_BesselTargets
        FILE ${QUICC_CURRENT_COMPONENT_LIB}_BesselTargets.cmake
        NAMESPACE ${QUICC_NAMESPACE}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/quicc
)
