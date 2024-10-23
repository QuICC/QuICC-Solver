# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Pointwise operators. This is done to avoid writing boiler
# plate code for each operator.

function(mapType2mlir TYPE MLIRTYPE)
    if(${TYPE} STREQUAL "double")
        set(${MLIRTYPE} "f64" PARENT_SCOPE)
    elseif(${TYPE} STREQUAL "std::complex<double>")
        set(${MLIRTYPE} "complexf64" PARENT_SCOPE)
    else()
        message(SEND_ERROR "unknown type")
    endif()
endfunction()

function(mapType2cuda TYPE CUDATYPE)
    if(${TYPE} STREQUAL "double")
        set(${CUDATYPE} "double" PARENT_SCOPE)
    elseif(${TYPE} STREQUAL "std::complex<double>")
        set(${CUDATYPE} "cuda::std::complex<double>" PARENT_SCOPE)
    else()
        message(SEND_ERROR "unknown type")
    endif()
endfunction()

function(camelCase STRIN STROUT)
    string (TOUPPER "${STRIN}" U)
    string (SUBSTRING "${U}"   0  1 A)
    string (SUBSTRING "${STRIN}" 1 -1 B)
    set (${STROUT} "${A}${B}" PARENT_SCOPE)
endfunction()

set(Ops "add;sub")
set(Types "double;std::complex<double>")
set(Layouts "DCCSC3D;S1CLCSC3D")

# Configure Shims
foreach(Op IN LISTS Ops)
    camelCase(${Op} OpCC)
    foreach(Layout IN LISTS Layouts)
        foreach(Type IN LISTS Types)
            if (Type STREQUAL "double" AND Layout STREQUAL "S1CLCSC3D")
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
