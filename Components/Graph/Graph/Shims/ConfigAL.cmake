# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC ALegendre operators. This is done to avoid writing boiler
# plate code for each operator.

# Here we follow the quiccir naming convention
# Ops: quiccir.jw.prj, quiccir.jw.int
set(Ops "prj;int")
# kinds: P, D1, ...
set(prjKinds "P;D1;Ll;LlD1;DivS1;DivS1Dp;LlDivS1;LlDivS1Dp")
set(intKinds "P;D1;Ll;LlD1;DivS1;DivS1Dp;LlDivS1;LlDivS1Dp;Ll2;DivLl;DivLlD1;DivLlDivS1;DivLlDivS1Dp")

# Special treatments
set(DivS1DpTreatment "diffPhi_m")
set(LlDivS1DpTreatment "diffPhi_m")
set(DivLlDivS1DpTreatment "diffPhi_m")

# Map op name to direction
set(prjOpDirection "bwd_t")
set(intOpDirection "fwd_t")

# Configure Shims
foreach(Op IN LISTS Ops)
    set(Kinds "${${Op}Kinds}")
    foreach(Kind IN LISTS Kinds)
        set(OpDirection "${${Op}OpDirection}")
        set(Treatment "${${Kind}Treatment}")
        if(NOT Treatment)
            set(Treatment "none_m")
        else()
            if("${OpDirection}" STREQUAL "bwd_t")
                set(Treatment "diffPhiPrj_m")
            else()
                set(Treatment "diffPhiInt_m")
            endif()
        endif()
        configure_file(
            "MlirALegendre${Op}Shims.cpp.in"
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/ALegendre/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
        target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
            PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/ALegendre/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
    endforeach()
endforeach()
