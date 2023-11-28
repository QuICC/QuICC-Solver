# Create component library
set(QUICC_CURRENT_SUBCOMPONENT_LIB Reduction)

# if cuda backend is present we compile as static lib
# to not propagate the dependency, otherwise header only lib
if(TARGET QuICC::Cuda)
  set(LIBTYPE "STATIC")
  set(INTTYPE "PUBLIC")
else()
  set(LIBTYPE "INTERFACE")
  set(INTTYPE "INTERFACE")
endif()

add_library(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} "${LIBTYPE}" "")

target_include_directories(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
    "${INTTYPE}"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

add_subdirectory(Cpu)
if(TARGET QuICC::Cuda)
  add_subdirectory(Cuda)
endif()

target_link_libraries(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  "${INTTYPE}"
    QuICC::View
    QuICC::Operator
    QuICC::Profiling
)


# Alias
add_library(${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${QUICC_CURRENT_SUBCOMPONENT_LIB}
  ALIAS ${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB})

# Export info
quicc_export_target(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  COMPONENT ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${QUICC_CURRENT_SUBCOMPONENT_LIB}
)

