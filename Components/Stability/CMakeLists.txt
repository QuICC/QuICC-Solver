# CMake setup for Stability component
set(QUICC_CURRENT_COMPONENT_DIR Components/Stability)
set(QUICC_CURRENT_COMPONENT_LIB Stability)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Components/Stability/cmake.d")

# Enable stability solver if SLEPc/PETSc are available
if(QUICC_HAVE_STABILITY_SOLVER)
  # Create component library
  add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
  set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)

  target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )

  set(PHS
    Stability/LinearStability.hpp
    Stability/StabilityBase.hpp
    Stability/MarginalCurve.hpp
    Stability/LinearStabilityFactory.hpp
  )

  # Add public header
  set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
    PROPERTIES
      PUBLIC_HEADER "${PHS}"
  )

  # Dependencies
  target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC
    QuICC::Io
    QuICC::Common
    QuICC::Framework
    QuICC::Profiling
    PkgConfig::SLEPC
    )

  # Add source directory
  add_subdirectory(Stability)

  # Add library alias (as the imported one)
  add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

  # Teststuite
  #include(cmake.d/setup/TestSuite.cmake)

  # Export info
  quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
    COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  )
endif()
