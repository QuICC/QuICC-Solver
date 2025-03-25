# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Cross operators. This is done to avoid writing boiler
# plate code for each operator.

set(Kinds "inertia;lorentz;induction")

# Configure Shims
foreach(Kind IN LISTS Kinds)
    configure_file(
        "MlirCrossShims.cpp.in"
        "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Cross/MlirShims/${Kind}.cpp"
    )
    target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
        PRIVATE
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Cross/MlirShims/${Kind}.cpp"
    )
endforeach()

