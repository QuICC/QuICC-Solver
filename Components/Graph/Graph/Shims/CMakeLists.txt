# Create component library
set(QUICC_CURRENT_SUBCOMPONENT_LIB Shims)
add_library(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} SHARED "")

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set(PHS
  ../../Graph/Shims/MlirShims.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    LINKER_LANGUAGE CXX
)

target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PRIVATE
    MlirShims.cpp
    MlirDotShims.cpp
    MlirTransposeShims.cpp
    MlirAllocShims.cpp
    MlirDeallocShims.cpp
)

if(TARGET QuICC::Cuda)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
    PRIVATE
      MlirAllocShims.cu
  )
endif()

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PUBLIC
    QuICC::ViewOps
)

target_link_options(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PUBLIC
    "-Wl,-no-as-needed"
)

# Generate Shims
include(ConfigFT.cmake)
include(ConfigAL.cmake)
include(ConfigJW.cmake)
include(ConfigPointwise.cmake)
include(ConfigMulConst.cmake)
include(ConfigCross.cmake)

# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${QUICC_CURRENT_SUBCOMPONENT_LIB}
  ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB})

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${QUICC_CURRENT_SUBCOMPONENT_LIB}
)

# Fix for clang and derived compilers
message(VERBOSE "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	target_compile_options(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} PRIVATE -fsized-deallocation)
endif()
