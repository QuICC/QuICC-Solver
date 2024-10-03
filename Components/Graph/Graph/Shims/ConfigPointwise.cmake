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

function(camelCase STRIN STROUT)
    string (TOUPPER "${STRIN}" U)
    string (SUBSTRING "${U}"   0  1 A)
    string (SUBSTRING "${STRIN}" 1 -1 B)
    set (${STROUT} "${A}${B}" PARENT_SCOPE)
endfunction()

set(Ops "add;sub")
set(Types "double;std::complex<double>")

# Configure Shims
foreach(Op IN LISTS Ops)
    camelCase(${Op} OpCC)
    foreach(Type IN LISTS Types)
        mapType2mlir(${Type} MlirType)
        configure_file(
            "MlirPointwiseShims.cpp.in"
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Pointwise/MlirShims/${Op}/${MlirType}${Backend}.cpp"
        )
        target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
            PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Pointwise/MlirShims/${Op}/${MlirType}${Backend}.cpp"
        )
    endforeach()
endforeach()
