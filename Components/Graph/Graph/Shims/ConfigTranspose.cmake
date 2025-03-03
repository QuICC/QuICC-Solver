# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Transpose operators. This is done to avoid writing boiler
# plate code for each operator.

include(ConfigUtils.cmake)

set(LayoutIns "DCCSC3D;")
set(LayoutOuts "DCCSC3D;")
set(TypeIns "std::complex<double>")
set(TypeOuts "std::complex<double>")

# Configure Shims
foreach(LayoutIn IN LISTS LayoutIns)
    foreach(LayoutOut IN LISTS LayoutOuts)
        foreach(TypeIn IN LISTS TypeIns)
            foreach(TypeOut IN LISTS TypeOuts)
                mapType2mlir(${TypeIn} MlirTypeIn)
                mapType2mlir(${TypeOut} MlirTypeOut)
                camelCase(${TypeIn} KindIn)
                camelCase(${TypeOut} KindOut)
                camelCase(${LayoutIn} LayoutIn)
                camelCase(${LayoutOut} LayoutOut)
                configure_file(
                    "MlirTransposeShims.cpp.in"
                    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${KindIn}${KindOut}${LayoutIn}${LayoutOut}.cpp"
                )
                target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
                    PRIVATE
                        "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${KindIn}${KindOut}${LayoutIn}${LayoutOut}.cpp"
                )
            endforeach()
        endforeach()
    endforeach()
endforeach()

