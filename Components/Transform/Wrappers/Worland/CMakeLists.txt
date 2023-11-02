set(Backends "Cpu")

# if(TARGET QuICC::Cuda)
#     list(APPEND Backends "Gpu")
# endif()

set(Kinds "Projector;Integrator")

set(ProjectorOps "P;D1;D1_P;D1R1;DivR1;DivR1_Zero;DivR1D1R1;DivR1D1R1_Zero;SphLapl;CylLaplh;CylLaplh_DivR1D1R1;D1CylLaplh;D1CylLaplh_D1DivR1D1R1;DivR1CylLaplh_Zero")
set(IntegratorOps "P;P_Zero;R1_Zero;DivR1_Zero;I2;I2_Zero;I2DivR1_Zero;I4DivR1_Zero;I6DivR1_Zero;DivR1D1R1_Zero;I2DivR1D1R1_Zero;I4DivR1D1R1_Zero")

set(IntegratorOpAttribute "CSL3D")
set(ProjectorOpAttribute "CSL3D")
set(IntegratorOpInOut "mods_t, phys_t")
set(ProjectorOpInOut "phys_t, mods_t")

# Ops
set(PPolyBuilder "Wnl")

set(P_ZeroPolyBuilder "Wnl")
set(P_ZeroMeanBuilder "void")

set(R1_ZeroPolyBuilder "rWnl")
set(R1_ZeroMeanBuilder "void")

set(D1PolyBuilder "dWnl")

set(D1_PPolyBuilder "dWnl")
set(D1_PMeanBuilder "Wnl")

set(D1R1PolyBuilder "drWnl")

set(ProjectorDivR1PolyBuilder "r_1Wnl<recurrence_t>")
set(IntegratorDivR1PolyBuilder "r_1Wnl<implicit_t>")

set(ProjectorDivR1_ZeroPolyBuilder "r_1Wnl<recurrence_t>")
set(IntegratorDivR1_ZeroPolyBuilder "r_1Wnl<implicit_t>")
set(DivR1_ZeroMeanBuilder "void")

set(ProjectorDivR1D1R1PolyBuilder "r_1drWnl<recurrence_t>")
set(IntegratorDivR1D1R1PolyBuilder "r_1drWnl<implicit_t>")

set(ProjectorDivR1D1R1_ZeroPolyBuilder "r_1drWnl<recurrence_t>")
set(IntegratorDivR1D1R1_ZeroPolyBuilder "r_1drWnl<implicit_t>")
set(DivR1D1R1_ZeroMeanBuilder "void")

set(SphLaplPolyBuilder "slaplWnl")

set(CylLaplhPolyBuilder "claplhWnl")

set(CylLaplh_DivR1D1R1PolyBuilder "claplhWnl")
set(CylLaplh_DivR1D1R1MeanBuilder "r_1drWnl<recurrence_t>")

set(D1CylLaplhPolyBuilder "dclaplhWnl")

set(D1CylLaplh_D1DivR1D1R1PolyBuilder "dclaplhWnl")
set(D1CylLaplh_D1DivR1D1R1MeanBuilder "dr_1drWnl")

set(DivR1CylLaplh_ZeroPolyBuilder "r_1claplhWnl")
set(DivR1CylLaplh_ZeroMeanBuilder "void")

set(I2PolyBuilder "Wnl")
set(I2INBuilder "I2")

set(I2_ZeroPolyBuilder "Wnl")
set(I2_ZeroINBuilder "I2")
set(I2_ZeroMeanBuilder "void")
set(I2_ZeroMeanINBuilder "void")

set(I2DivR1_ZeroPolyBuilder "r_1Wnl<implicit_t>")
set(I2DivR1_ZeroINBuilder "I2")
set(I2DivR1_ZeroMeanBuilder "void")
set(I2DivR1_ZeroMeanINBuilder "void")

set(I2DivR1D1R1_ZeroPolyBuilder "r_1drWnl<implicit_t>")
set(I2DivR1D1R1_ZeroINBuilder "I2")
set(I2DivR1D1R1_ZeroMeanBuilder "void")
set(I2DivR1D1R1_ZeroMeanINBuilder "void")

set(I4DivR1_ZeroPolyBuilder "r_1Wnl<implicit_t>")
set(I4DivR1_ZeroINBuilder "I4")
set(I4DivR1_ZeroMeanBuilder "void")
set(I4DivR1_ZeroMeanINBuilder "void")

set(I4DivR1D1R1_ZeroPolyBuilder "r_1drWnl<implicit_t>")
set(I4DivR1D1R1_ZeroINBuilder "I4")
set(I4DivR1D1R1_ZeroMeanBuilder "void")
set(I4DivR1D1R1_ZeroMeanINBuilder "void")

set(I6DivR1_ZeroPolyBuilder "r_1Wnl<implicit_t>")
set(I6DivR1_ZeroINBuilder "I6")
set(I6DivR1_ZeroMeanBuilder "void")
set(I6DivR1_ZeroMeanINBuilder "void")

# cpu/gpu
set(CpuMemRes "Cpu::NewDelete")
set(GpuMemRes "Cuda::Malloc")

set(GpuLayout "#define QUICC_JW_ROW_MAJOR")

set(GpuHppHeaders "#include <cuda/std/complex>
#include \"Memory/Memory.hpp\""
)

set(GpuCppHeaders "#include \"Memory/Cuda/Malloc.hpp\"
#include \"Cuda/CudaUtil.hpp\""
)

set(CpuBackendNamespace "Cpu")
set(GpuBackendNamespace "Cuda")

set(GpuStorage
"/// @brief Device data pointers
    mutable QuICC::Memory::MemBlock<std::uint32_t> mDataPointers_d;

    /// @brief Device data indices
    mutable QuICC::Memory::MemBlock<std::uint32_t> mDataIndices_d;

    /// @brief Temporary storage for flattened output on device
    mutable QuICC::Memory::MemBlock<cuda::std::complex<double>> mModsFlat_d;

    /// @brief Temporary storage for flattened input on device
    mutable QuICC::Memory::MemBlock<cuda::std::complex<double>> mPhysFlat_d;"
)

set(CpuSetDataViews
"// set pointers
    ViewBase<std::uint32_t> dataPointers[rank];
    dataPointers[1] = ViewBase<std::uint32_t>(mDataPointers.data(), mDataPointers.size());
    ViewBase<std::uint32_t> dataIndices[rank];
    dataIndices[1] = ViewBase<std::uint32_t>(mDataIndices.data(), mDataIndices.size());

    // set views
    mModsView = mods_t(mModsFlat.data(), mModsFlat.size(), modsDims.data(), dataPointers, dataIndices);
    mPhysView = phys_t(mPhysFlat.data(), mPhysFlat.size(), physDims.data(), dataPointers, dataIndices);"
)

set(GpuSetDataViews
"// alloc device memory
    mModsFlat_d = std::move(QuICC::Memory::MemBlock<cuda::std::complex<double>>(mModsFlat.size(), mem.get()));
    mDataPointers_d = std::move(QuICC::Memory::MemBlock<std::uint32_t>(mDataPointers.size(), mem.get()));
    mDataIndices_d = std::move(QuICC::Memory::MemBlock<std::uint32_t>(mDataIndices.size(), mem.get()));
    mPhysFlat_d = std::move(QuICC::Memory::MemBlock<cuda::std::complex<double>>(mPhysFlat.size(), mem.get()));

    // set device pointers
    ViewBase<std::uint32_t> dataPointers_d[rank];
    dataPointers_d[1] = ViewBase<std::uint32_t>(mDataPointers_d.data(), mDataPointers_d.size());
    ViewBase<std::uint32_t> dataIndices_d[rank];
    dataIndices_d[1] = ViewBase<std::uint32_t>(mDataIndices_d.data(), mDataIndices_d.size());

    // copy to device memory
    cudaErrChk(cudaMemcpyAsync(reinterpret_cast<std::complex<double>*>(mModsFlat_d.data()), mModsFlat.data(),
        mModsFlat_d.size()*sizeof(cuda::std::complex<double>), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpyAsync(mDataPointers_d.data(), mDataPointers.data(),
        mDataPointers_d.size()*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpyAsync(mDataIndices_d.data(), mDataIndices.data(),
        mDataIndices_d.size()*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpyAsync(reinterpret_cast<std::complex<double>*>(mPhysFlat_d.data()), mPhysFlat.data(),
        mPhysFlat_d.size()*sizeof(cuda::std::complex<double>), cudaMemcpyHostToDevice));

    // set device views
    mModsView = mods_t(reinterpret_cast<std::complex<double>*>(mModsFlat_d.data()), mModsFlat_d.size(), modsDims.data(), dataPointers_d, dataIndices_d);
    mPhysView = phys_t(reinterpret_cast<std::complex<double>*>(mPhysFlat_d.data()), mPhysFlat_d.size(), physDims.data(), dataPointers_d, dataIndices_d);"
)

set(GpuIntegratorHostToDevice
"cudaErrChk(cudaMemcpy(mPhysFlat_d.data(), mPhysFlat.data(),
    mPhysFlat_d.size()*sizeof(cuda::std::complex<double>), cudaMemcpyHostToDevice));"
)

set(GpuIntegratorDeviceToHost
"cudaErrChk(cudaMemcpy(mModsFlat.data(), mModsFlat_d.data(),
    mModsFlat.size()*sizeof(cuda::std::complex<double>), cudaMemcpyDeviceToHost));"
)

set(GpuProjectorHostToDevice
"cudaErrChk(cudaMemcpy(mModsFlat_d.data(), mModsFlat.data(),
    mModsFlat_d.size()*sizeof(cuda::std::complex<double>), cudaMemcpyHostToDevice));"
)

set(GpuProjectorDeviceToHost
"cudaErrChk(cudaMemcpy(mPhysFlat.data(), mPhysFlat_d.data(),
    mPhysFlat.size()*sizeof(cuda::std::complex<double>), cudaMemcpyDeviceToHost));"
)

# end setup

foreach(Backend IN LISTS Backends)
    set(Tag "view${Backend}_t")
    set(Layout "${${Backend}Layout}")
    set(HppHeaders "${${Backend}HppHeaders}")
    set(Storage "${${Backend}Storage}")
    set(CppHeaders "${${Backend}CppHeaders}")
    set(BackendNamespace "${${Backend}BackendNamespace}")
    set(MemRes "${${Backend}MemRes}")
    set(SetDataViews "${${Backend}SetDataViews}")
    foreach(Kind IN LISTS Kinds)
        set(Ops "${${Kind}Ops}")
        set(OpAttribute "${${Kind}OpAttribute}")
        set(OpInOut "${${Kind}OpInOut}")
        set(HostToDevice "${${Backend}${Kind}HostToDevice}")
        set(DeviceToHost "${${Backend}${Kind}DeviceToHost}")
        foreach(Op IN LISTS Ops)
            set(PolyBuilder "${${Op}PolyBuilder}")
            if (NOT PolyBuilder)
                set(PolyBuilder "${${Kind}${Op}PolyBuilder}")
            endif()

            set(OpBuilder "Operator<${PolyBuilder}>")
            set(MeanBuilder "${${Op}MeanBuilder}")
            set(INBuilder "${${Op}INBuilder}")
            if(MeanBuilder)
                if (NOT (MeanBuilder STREQUAL "void"))
                    set(MeanHeader "#include \"QuICC/Polynomial/Worland/${MeanBuilder}.hpp\"")
                endif()
                set(OpBuilder "OperatorWithMean<${PolyBuilder}, ${MeanBuilder}>")
            endif()
            if(INBuilder)
                set(INHeader "#include \"QuICC/SparseSM/Worland/${INBuilder}.hpp\"")
                set(OpBuilder "OperatorIN<${PolyBuilder}, ::QuICC::SparseSM::Worland::${INBuilder}>")
            endif()
            if(INBuilder AND MeanBuilder)
                set(INHeader "#include \"QuICC/SparseSM/Worland/${INBuilder}.hpp\"")
                set(OpBuilder "OperatorINWithMean<${PolyBuilder}, ::QuICC::SparseSM::Worland::${INBuilder}>")
            endif()


            set(PolyBuilderScale "${${Op}PolyBuilderScale}")

            configure_file(
                "OpV.hpp.in"
                "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Worland/${Kind}/${Op}${Tag}.hpp.inc"
            )
            set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
                PROPERTIES
                    PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Worland/${Kind}/${Op}${Tag}.hpp.inc"
            )

            configure_file(
                "${Kind}/OpV.cpp.in"
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Worland/${Kind}/${Op}${Tag}.cpp"
            )
            target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Worland/${Kind}/${Op}${Tag}.cpp"
            )

        endforeach()
    endforeach()
endforeach()
