# CMake setup for Common component
set(QUICC_CURRENT_COMPONENT_DIR Components/Common)
set(QUICC_CURRENT_COMPONENT_LIB ${QUICC_COMMON_COMPONENT_LIB})

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

include(cmake.d/setup/RegisterArithmetics.cmake)
include(cmake.d/setup/RegisterModelOperator.cmake)
include(cmake.d/setup/RegisterModelOperatorBoundary.cmake)
include(cmake.d/setup/RegisterNonDimensional.cmake)
include(cmake.d/setup/RegisterPhysicalNames.cmake)
include(cmake.d/setup/RegisterPseudospectralTag.cmake)
include(cmake.d/setup/RegisterRuntimeStatus.cmake)
include(cmake.d/setup/RegisterSolveTiming.cmake)

# Create component library
add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"  )

target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" QuICC::Profiling::Timers)

# Add MPI dependency
if(QUICC_MPI)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" MPI::MPI_CXX)
endif(QUICC_MPI)

# Boost
if(QUICC_MULTPRECISION)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" Boost::headers)
endif()

# Eigen
if(TARGET Eigen3::Eigen)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" Eigen3::Eigen)
endif()

# UMFPACK
if(TARGET UMFPACK::UMFPACK)
  target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" UMFPACK::UMFPACK)
endif()

target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" hdf5::hdf5)

# Add source directory
add_subdirectory(src)

# Teststuite
include(cmake.d/setup/TestSuite.cmake)

# Export
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT Common
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC
    )
