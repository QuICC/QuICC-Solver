# Note, new ops don't necessarily need to be generated

# Here we follow the quiccir naming convention
# Ops: quiccir.jw.prj, quiccir.jw.int
set(Ops "prj;int")
# kinds: P, D1, ...
set(prjKinds "P;D1")

# QuICC Op Builders
set(prjOpDirection "bwd_t")
set(PPolyBuilder "Wnl")
set(D1PolyBuilder "dWnl")

# Configure Shims
foreach(Op IN LISTS Ops)
    set(Kinds "${${Op}Kinds}")
    foreach(Kind IN LISTS Kinds)
        set(OpDirection "${${Op}OpDirection}")
        set(PolyBuilder "${${Kind}PolyBuilder}")
        set(OpBuilder "Operator<${PolyBuilder}>")

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