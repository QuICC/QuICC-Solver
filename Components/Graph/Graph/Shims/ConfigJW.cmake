# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Worland operators. This is done to avoid writing boiler
# plate code for each operator.

# Here we follow the quiccir naming convention
# Ops: quiccir.jw.prj, quiccir.jw.int
set(Ops "prj;int")
# kinds: P, D1, ...
set(prjKinds "P;D1;D1R1;D1_P;DivR1;DivR1_Zero;DivR1D1R1;DivR1D1R1_Zero;SphLapl;CylLaplh;CylLaplh_DivR1D1R1;D1CylLaplh;D1CylLaplh_D1DivR1D1R1;DivR1CylLaplh_Zero;Energy;EnergyR2;EnergyD1R1;EnergySLaplR2;RadialPower;RadialPowerDivR1;RadialPowerDivR1D1R1")
set(intKinds "P;P_Zero;R1_Zero;DivR1_Zero;I2;I2_Zero;I2DivR1_Zero;I4DivR1_Zero;I6DivR1_Zero;DivR1D1R1_Zero;I2DivR1D1R1_Zero;I4DivR1D1R1_Zero;Energy;EnergyR2;EnergyD1R1;EnergySLaplR2;RadialPower;RadialPowerDivR1;RadialPowerDivR1D1R1")
# map op name to direction
set(prjOpDirection "bwd_t")
set(intOpDirection "fwd_t")

# Configure Shims
foreach(Op IN LISTS Ops)
    set(Kinds "${${Op}Kinds}")
    foreach(Kind IN LISTS Kinds)
        set(OpDirection "${${Op}OpDirection}")
        configure_file(
            "MlirWorland${Op}Shims.cpp.in"
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Worland/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
        target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
            PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Worland/MlirShims/${Op}/${Kind}${Backend}.cpp"
        )
    endforeach()
endforeach()
