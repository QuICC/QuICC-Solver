set(Backends "Cpu")

if(TARGET QuICC::Cuda)
    list(APPEND Backends "Gpu")
endif()

set(Kinds "Projector;Integrator")

set(ProjectorOps "P;D1;Ll;LlD1;DivS1;DivS1Dp;LlDivS1;LlDivS1Dp")
set(IntegratorOps "P;D1;Ll;LlD1;DivS1;DivS1Dp;LlDivS1;LlDivS1Dp;Ll2;DivLl;DivLlD1;DivLlDivS1;DivLlDivS1Dp")

set(IntegratorOpAttribute "S1CLCSC3D")
set(ProjectorOpAttribute "CS1RL3D")
set(IntegratorOpInOut "mods_t, phys_t")
set(ProjectorOpInOut "phys_t, mods_t")

# Ops
set(DivS1DpTreatment "diffPhi_m")
set(LlDivS1DpTreatment "diffPhi_m")
set(DivLlDivS1DpTreatment "diffPhi_m")

# cpu/gpu
set(CpuMemRes "Cpu::NewDelete")
set(GpuMemRes "Cuda::Malloc")

set(GpuLayout "#define QUICC_AL_ALL_ROW_MAJOR")

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
        if("${Kind}" STREQUAL "Projector")
            set(OpDirection "bwd_t")
        else()
            set(OpDirection "fwd_t")
        endif()
        foreach(Op IN LISTS Ops)
            set(PolyBuilder "${${Op}PolyBuilder}")
            set(Diff "${${Op}Diff}")
            set(PolyBuilderScale "${${Op}PolyBuilderScale}")
            set(Treatment "${${Op}Treatment}")
            if(NOT Treatment)
                set(Treatment "none_m")
            else()
                if("${Kind}" STREQUAL "Projector")
                    set(Treatment "diffPhiPrj_m")
                else()
                    set(Treatment "diffPhiInt_m")
                endif()
            endif()

            configure_file(
                "OpV.hpp.in"
                "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/ALegendre/${Kind}/${Op}${Tag}.hpp.inc"
            )
            set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
                PROPERTIES
                    PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/ALegendre/${Kind}/${Op}${Tag}.hpp.inc"
            )

            configure_file(
                "${Kind}/OpV.cpp.in"
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/ALegendre/${Kind}/${Op}${Tag}.cpp"
            )
            target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/ALegendre/${Kind}/${Op}${Tag}.cpp"
            )

        endforeach()
    endforeach()
endforeach()
