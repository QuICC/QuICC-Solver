# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Worland operators. This is done to avoid writing boiler
# plate code for each operator.

# Here we follow the quiccir naming convention
# Ops: quiccir.jw.prj, quiccir.jw.int
set(Ops "prj;int")
# kinds: P, D1, ...
set(prjKinds "P;D1;D1R1")

# map op name to direction
set(prjOpDirection "bwd_t")
set(intOpDirection "fwd_t")

# Configure Shims
foreach(Op IN LISTS Ops)
    set(Kinds "${${Op}Kinds}")
    foreach(Kind IN LISTS Kinds)
        set(OpDirection "${${Op}OpDirection}")
        configure_file(
            "MlirWorlandShims.cpp.in"
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Worland/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
        target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
            PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Worland/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
    endforeach()
endforeach()