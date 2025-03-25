# This scripts generates the Shims (glue) between the mlir library call
# and the QuICC Mul Const operators. This is done to avoid writing boiler
# plate code for each operator.

set(Ops "Theta;R")
set(ThetaKinds "buoyancy;transport")
set(ThetaGridBuilder "QuICC::Polynomial::Quadrature::WorlandRule")
set(ThetaDir "2")
set(RKinds "coriolis_cos;coriolis_sin")
set(RGridBuilder "QuICC::Polynomial::Quadrature::LegendreRule")
set(RDir "1")


# Configure Shims
foreach(Op IN LISTS Ops)
    set(Kinds "${${Op}Kinds}")
    foreach(Kind IN LISTS Kinds)
        set(Dir "${${Op}Dir}")
        set(GridBuilder "${${Op}GridBuilder}")
        set(Functor "MulRFunctor")
        if (Kind STREQUAL "coriolis_sin")
            set(Functor "MulSinFunctor")
        endif()
        if (Kind STREQUAL "coriolis_cos")
            set(Functor "MulCosFunctor")
        endif()
        configure_file(
            "MlirMulConstShims.cpp.in"
            "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Slicewise/MlirShims/${Op}/${Kind}.cpp"
        )
        target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
            PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/Slicewise/MlirShims/${Op}/${Kind}.cpp"
        )
    endforeach()
endforeach()
