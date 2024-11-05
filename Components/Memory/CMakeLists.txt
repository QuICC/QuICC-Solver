# CMake setup for Memory utils
set(QUICC_CURRENT_COMPONENT_LIB Memory)
set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QUICC_CURRENT_COMPONENT_LIB}::")
message(DEBUG "QUICC_CURRENT_COMPONENT_LIB_NAMESPACE: ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}" )


add_library(${QUICC_CURRENT_COMPONENT_LIB} STATIC
  Memory/Cpu/NewDelete.cpp
  Memory/Cpu/Pool.cpp
  Memory/Cuda/Malloc.cu
)

set_property(TARGET ${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# set(PHS
#   Memory/Memory.hpp
#   Memory/MemoryResource.hpp
#   Memory/Pensieve.hpp
#   Memory/Cpu/NewDelete.hpp
#   Memory/Cpu/Pool.hpp
#   Memory/Cuda/Malloc.hpp
# )

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

# Fix for clang and derived compilers
message(VERBOSE "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	target_compile_options(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE -fsized-deallocation)
endif()


# Cuda
if(TARGET QuICC::Cuda)
    target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB} QuICC::Cuda)
endif()


# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB} ALIAS "${QUICC_CURRENT_COMPONENT_LIB}")

# Teststuite
option(QUICC_TESTSUITE_MEMORY "Enable Memory component testsuite?" OFF)
if(QUICC_TESTSUITE_MEMORY)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_MEMORY)

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB}
  DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/Memory
    )
