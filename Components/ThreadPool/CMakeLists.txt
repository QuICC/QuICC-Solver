option(QUICC_USE_THREADPOOL "Enable threadpool component?" OFF)
set(QUICC_USE_THREADPOOL ${QUICC_USE_THREADPOOL} PARENT_SCOPE)
if(NOT QUICC_USE_THREADPOOL)
  return()
endif()

# CMake setup for ThreadPool
set(QUICC_CURRENT_COMPONENT_LIB ThreadPool)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )

# Set library sources visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

add_library(${QUICC_CURRENT_COMPONENT_LIB} "")
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB} PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    QUICC_USE_THREADPOOL
)

set(PHS
  ThreadPool/ThreadPool.hpp
  ThreadPool/ThreadPoolBoost.hpp
  ThreadPool/QuICCThreads.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Boost with thread component
include(QuICCFindBoostRequired)
quicc_find_boost_required(COMPONENTS thread)
target_link_libraries("${QUICC_CURRENT_COMPONENT_LIB}" PUBLIC Boost::headers Boost::thread)

add_subdirectory(ThreadPool)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS
  "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
option(QUICC_TESTSUITE_THREADPOOL "Enable ThreadPool component testsuite?" OFF)
if(QUICC_TESTSUITE_THREADPOOL)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_THREADPOOL)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
