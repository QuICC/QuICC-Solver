# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Transpose operators. This is done to avoid writing boiler
# plate code for each operator.

include(ConfigUtils.cmake)

set(LayoutIns "DCCSC3D")
set(LayoutOuts "DCCSC3D")
set(TypeIns "std::complex<double>")
set(TypeOuts "std::complex<double>")

set(GroupSizes "1;2")

# Configure Shims
foreach(LayoutIn IN LISTS LayoutIns)
    foreach(LayoutOut IN LISTS LayoutOuts)
        foreach(TypeIn IN LISTS TypeIns)
            foreach(TypeOut IN LISTS TypeOuts)
                mapType2mlir(${TypeIn} MlirTypeIn)
                mapType2mlir(${TypeOut} MlirTypeOut)

                foreach(GroupSize IN LISTS GroupSizes)
                    math(EXPR GroupSizeM1 "${GroupSize}-1")
                    set(FunName "_ciface_quiccir_transpose")
                    foreach(It RANGE 0 ${GroupSizeM1})
                        string(APPEND FunName "_${MlirTypeOut}_${LayoutOut}")
                    endforeach()
                    foreach(It RANGE 0 ${GroupSizeM1})
                        string(APPEND FunName "_${MlirTypeIn}_${LayoutIn}")
                    endforeach()

                    set(FunSignature "${FunName}(void* obj")

                    foreach(It RANGE 0 ${GroupSizeM1})
                        string(APPEND FunSignature ", ViewDescriptor<${TypeOut}, std::uint32_t, 3>* pOut${It}")
                    endforeach()

                    foreach(It RANGE 0 ${GroupSizeM1})
                        string(APPEND FunSignature ", const ViewDescriptor<${TypeIn}, std::uint32_t, 3>* pIn${It}")
                    endforeach()
                    string(APPEND FunSignature ")")

                    # View defs

                    configure_file(
                        "MlirTransposeShims.cpp.in"
                        "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${FunName}.cpp"
                    )
                    target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
                        PRIVATE
                            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${FunName}.cpp"
                    )
                endforeach(GroupSize IN LISTS GroupSizes)


            endforeach()
        endforeach()
    endforeach()
endforeach()

