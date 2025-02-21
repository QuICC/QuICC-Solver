message(STATUS "... adding stability tools")

# CMake setup for Stability tool
set(QUICC_CURRENT_TOOL Stability)
set(QUICC_CURRENT_TOOL_DIR Tools/${QUICC_CURRENT_TOOL})
set(QUICC_CURRENT_TOOL_LIB Tools_${QUICC_CURRENT_TOOL})
set(QUICC_CURRENT_TOOL_LIB_NAMESPACE "Tools::${QUICC_CURRENT_TOOL}::")

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${QUICC_CURRENT_TOOL_DIR}/cmake.d")

#  Generate specific nondimension parameters
include(cmake.d/setup/RegisterNonDimensional.cmake)

# Enable stability solver if SLEPc/PETSc are available
if(QUICC_HAVE_STABILITY_SOLVER)
  message(STATUS "SETTING UP STABILITY TOOL")
  # Create tool library
  add_library(${QUICC_CURRENT_TOOL_LIB} "")
  set_target_properties(${QUICC_CURRENT_TOOL_LIB} PROPERTIES LINKER_LANGUAGE CXX)

  target_include_directories(${QUICC_CURRENT_TOOL_LIB}
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
  set_target_properties(${QUICC_CURRENT_TOOL_LIB}
    PROPERTIES
      PUBLIC_HEADER "${PHS}"
  )

  # Dependencies
  target_link_libraries(${QUICC_CURRENT_TOOL_LIB} PUBLIC
    QuICC::Io
    QuICC::Common
    QuICC::Framework
    QuICC::Profiling
    PkgConfig::SLEPC
    )

  # Add source directory
  add_subdirectory(Stability)

  # Add library alias (as the imported one)
  add_library(${QUICC_NAMESPACE}Tools::${QUICC_CURRENT_TOOL}
    ALIAS "${QUICC_CURRENT_TOOL_LIB}")

  # Teststuite
  #include(cmake.d/setup/TestSuite.cmake)

  # Export info
  quicc_export_target(${QUICC_CURRENT_TOOL_LIB}
    COMPONENT ${QUICC_CURRENT_TOOL_LIB}
    DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/Stability"
  )
endif()
