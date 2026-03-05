# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Transpose operators. This is done to avoid writing boiler
# plate code for each operator.

include(ConfigUtils.cmake)

set(Backends "Cpu")
if(TARGET QuICC::Cuda)
if(QUICC_USE_PFSOLVE)
else()
    list(APPEND Backends "Cuda")
endif()
endif()

# These 3 lists will be iterated over together
set(LayoutOutsCpu "DCCSC3D;DCCSC3D;DCCSC3D;S1CLCSC3D")
set(LayoutInsCpu "DCCSC3D;DCCSC3D;S1CLCSC3D;DCCSC3D")
if(QUICC_USE_PFSOLVE)
#set(LayoutOutsCuda "DCCSC3D;DCCSC3D;DCCSC3D;S1CLCSC3D")
#set(LayoutInsCuda "DCCSC3D;DCCSC3D;S1CLCSC3D;DCCSC3D")
#set(LayoutOutsCuda "DCCSC3DJIK;DCCSC3D;DCCSC3DJIK;S1CLCSC3DJIK")
#set(LayoutInsCuda "DCCSC3D;DCCSC3DJIK;S1CLCSC3DJIK;DCCSC3DJIK")
else()
set(LayoutOutsCuda "DCCSC3DJIK;DCCSC3D;DCCSC3DJIK;S1CLCSC3DJIK")
set(LayoutInsCuda "DCCSC3D;DCCSC3DJIK;S1CLCSC3DJIK;DCCSC3DJIK")
endif()
set(Perms "201;120;201;120")

# These 3 lists will be iterated over independently
set(TypeIns "std::complex<double>")
set(TypeOuts "std::complex<double>")
set(GroupSizes "1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16")

# Configure Shims
foreach(Backend IN LISTS Backends)
    set(LayoutOuts "${LayoutOuts${Backend}}")
    set(LayoutIns "${LayoutIns${Backend}}")
    list(LENGTH Perms len)
    math(EXPR lenM1 "${len} - 1")
    foreach(ItP RANGE 0 ${lenM1})
        list(GET Perms ${ItP} Perm)
        list(GET LayoutIns ${ItP} LayoutIn)
        list(GET LayoutOuts ${ItP} LayoutOut)
        foreach(TypeIn IN LISTS TypeIns)
            foreach(TypeOut IN LISTS TypeOuts)
                mapType2mlir(${TypeIn} MlirTypeIn)
                mapType2mlir(${TypeOut} MlirTypeOut)

                foreach(GroupSize IN LISTS GroupSizes)
                    math(EXPR GroupSizeM1 "${GroupSize}-1")
                    set(FunName "_ciface_quiccir_transpose_${Perm}")
                    set(FileName "${FunName}_${GroupSize}_${LayoutOut}_${LayoutIn}")
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

                    # Set lds
                    # needef for padding of fourier ops
                    if(LayoutIn STRGREATER_EQUAL "DCCSC3D" AND LayoutOut STRGREATER_EQUAL "DCCSC3D")
                        if (Perm STREQUAL "201")
                            set(Lds "std::uint32_t lds = pIn0->dataSize / pIn0->cooSize;")
                        else()
                            set(Lds "std::uint32_t lds = pOut0->dataSize / pOut0->cooSize;")
                        endif()
                    else()
                        set(Lds "")
                    endif()


                    # View defs
                    set(ViewDefs "")
                    foreach(It RANGE 0 ${GroupSizeM1})
                        if(LayoutIn STRGREATER_EQUAL "DCCSC3D" AND
                        LayoutOut STRGREATER_EQUAL "DCCSC3D" AND
                        Perm STREQUAL "201")
                            string(APPEND ViewDefs "    Tin::value_type viewIn${It}(pIn${It}->data, pIn${It}->dataSize, pIn${It}->dims, pointersIn, indicesIn, lds);\n")
                        else()
                        string(APPEND ViewDefs "    Tin::value_type viewIn${It}(pIn${It}->data, pIn${It}->dataSize, pIn${It}->dims, pointersIn, indicesIn);\n")
                        endif()
                    endforeach()
                    foreach(It RANGE 0 ${GroupSizeM1})
                        if(LayoutIn STRGREATER_EQUAL "DCCSC3D" AND
                        LayoutOut STRGREATER_EQUAL "DCCSC3D" AND
                        Perm STREQUAL "120")
                            string(APPEND ViewDefs "    Tout::value_type viewOut${It}(pOut${It}->data, pOut${It}->dataSize, pOut${It}->dims, pointersOut, indicesOut, lds);\n")
                        else()
                            string(APPEND ViewDefs "    Tout::value_type viewOut${It}(pOut${It}->data, pOut${It}->dataSize, pOut${It}->dims, pointersOut, indicesOut);\n")
                        endif()
                    endforeach()

                    # VecView defs
                    set(VecViewDefs "    std::vector<Tin::value_type> viewIns = {")
                    foreach(It RANGE 0 ${GroupSizeM1})
                        if(It EQUAL GroupSizeM1)
                            string(APPEND VecViewDefs "viewIn${It}")
                        else()
                            string(APPEND VecViewDefs "viewIn${It}, ")
                        endif()
                    endforeach()
                    string(APPEND VecViewDefs "};\n")
                    string(APPEND VecViewDefs "    std::vector<Tout::value_type> viewOuts = {")
                    foreach(It RANGE 0 ${GroupSizeM1})
                        if(It EQUAL GroupSizeM1)
                            string(APPEND VecViewDefs "viewOut${It}")
                        else()
                            string(APPEND VecViewDefs "viewOut${It}, ")
                        endif()
                    endforeach()
                    string(APPEND VecViewDefs "};\n")

                    # Configure file
                    configure_file(
                        "MlirTransposeShims.cpp.in"
                        "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${Backend}/${FileName}.cpp"
                    )
                    target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
                        PRIVATE
                            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Transpose/MlirShims/${Backend}/${FileName}.cpp"
                    )
                endforeach(GroupSize IN LISTS GroupSizes)
            endforeach()
        endforeach()
    endforeach()
endforeach(Backend IN LISTS Backends)



