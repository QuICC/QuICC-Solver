# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Pointwise operators. This is done to avoid writing boiler
# plate code for each operator.

include(ConfigUtils.cmake)

set(Ops "add;sub")
set(Types "double;std::complex<double>")
set(Layouts "DCCSC3D;DCCSC3DJIK;S1CLCSC3D;S1CLCSC3DJIK")

# Configure Shims
foreach(Op IN LISTS Ops)
    camelCase(${Op} OpCC)
    foreach(Layout IN LISTS Layouts)
        foreach(Type IN LISTS Types)
            # These layouts appear only in spectral space
            if (Type STREQUAL "double" AND
                (Layout STREQUAL "S1CLCSC3D"
                    OR Layout STREQUAL "S1CLCSC3DJIK"
                    OR Layout STREQUAL "DCCSC3DJIK"))
                continue()
            endif()
            mapType2mlir(${Type} MlirType)
            mapType2cuda(${Type} CudaType)
            configure_file(
                "MlirPointwiseShims.cpp.in"
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Pointwise/MlirShims/${Op}/${MlirType}${Layout}${Backend}.cpp"
            )
            target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
                PRIVATE
                    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Pointwise/MlirShims/${Op}/${MlirType}${Layout}${Backend}.cpp"
            )
        endforeach()
    endforeach()
endforeach()
