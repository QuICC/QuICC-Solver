set(Backends "Cpu")

if(TARGET QuICC::Cuda)
    list(APPEND Backends "Gpu")
endif()

set(Kinds "Energy")

set(EnergyOps "Energy;EnergyR2;EnergyD1R1;EnergySLaplR2")

set(EnergyOpAttribute "CSL3D")

set(EnergySizeRow "this->mspSetup->blockSize();")
set(EnergySizeCol "1;")

# cpu/gpu
set(CpuMemRes "Cpu::NewDelete")
set(GpuMemRes "Cuda::Malloc")

set(GpuLayout "#define QUICC_FD_SPHERE_ALL_ROW_MAJOR")

set(GpuHppHeaders
"
#include <cuda/std/complex>
#include \"Memory/Memory.hpp\"
"
)

set(GpuCppHeaders
"
#include \"Memory/Cuda/Malloc.hpp\"
#include \"Cuda/CudaUtil.hpp\"
#include \"ViewOps/Pointwise/Cuda/Pointwise.hpp\"
#include \"ViewOps/Reduction/Cuda/Reduction.hpp\"
"
)

set(CpuBackendNamespace "Cpu")
set(GpuBackendNamespace "Cuda")

set(GpuStorage
"
    /// @brief Device data pointers
    mutable QuICC::Memory::MemBlock<std::uint32_t> mDataPointers_d;

    /// @brief Device data indices
    mutable QuICC::Memory::MemBlock<std::uint32_t> mDataIndices_d;

    /// @brief Temporary storage for flattened output on device
    mutable QuICC::Memory::MemBlock<cuda::std::complex<double>> mModsFlat_d;

    /// @brief Temporary storage for flattened temp/output on device
    mutable QuICC::Memory::MemBlock<double> mModsAbsFlat_d;

    /// @brief Temporary storage for flattened temp/output on device
    mutable QuICC::Memory::MemBlock<double> mModsRedFlat_d;

    /// @brief Temporary storage for flattened input on device
    mutable QuICC::Memory::MemBlock<cuda::std::complex<double>> mPhysFlat_d;

    /// @brief Temporary storage for flattened temp/output on device
    mutable QuICC::Memory::MemBlock<double> mPhysAbsFlat_d;
"
)

set(CpuSetDataViews
"// set pointers
    ViewBase<std::uint32_t> dataPointers[rank];
    dataPointers[1] = ViewBase<std::uint32_t>(mDataPointers.data(), mDataPointers.size());
    ViewBase<std::uint32_t> dataIndices[rank];
    dataIndices[1] = ViewBase<std::uint32_t>(mDataIndices.data(), mDataIndices.size());

    // set views
    mModsView = mods_t(mModsFlat.data(), mModsFlat.size(), modsDims.data(), dataPointers, dataIndices);
    mPhysView = phys_t(mPhysFlat.data(), mPhysFlat.size(), physDims.data(), dataPointers, dataIndices);

    mPhysAbsFlat.resize(M * nCols);
    mPhysAbsView = physAbs_t(mPhysAbsFlat.data(), mPhysAbsFlat.size(), physDims.data(), dataPointers, dataIndices);

    std::array<std::uint32_t, rank> redDims = {1, N, P};
    mModsRedFlat.resize(nCols);
    mModsRedView = modsRed_t(mModsRedFlat.data(), mModsRedFlat.size(), redDims.data(), dataPointers, dataIndices);"
)

set(GpuSetDataViews
"
    // alloc device memory
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
    mPhysView = phys_t(reinterpret_cast<std::complex<double>*>(mPhysFlat_d.data()), mPhysFlat_d.size(), physDims.data(), dataPointers_d, dataIndices_d);
"
)

set(PowerGpuSetPowerView
"
    mModsAbsFlat.resize(inFlatSize);

    // alloc device mem
    mModsAbsFlat_d = std::move(QuICC::Memory::MemBlock<double>(mModsAbsFlat.size(), mem.get()));

    // set device views
    mModsAbsView = modsAbs_t(mModsAbsFlat_d.data(), mModsAbsFlat_d.size(), modsDims.data(), dataPointers_d, dataIndices_d);
")

set(EnergyGpuSetPowerView
"
    // alloc device mem
    mModsAbsFlat_d = std::move(QuICC::Memory::MemBlock<double>(inFlatSize, mem.get()));

    // set device views
    mModsAbsView = modsAbs_t(mModsAbsFlat_d.data(), mModsAbsFlat_d.size(), modsDims.data(), dataPointers_d, dataIndices_d);

    // setup reductor
    mRed = std::make_unique<QuICC::Reduction::Cuda::Op<modsRed_t, modsAbs_t, 0>>(mem);

    // set pointers
    ViewBase<std::uint32_t> dataPointers2D[rank-1];
    dataPointers2D[0] = ViewBase<std::uint32_t>(mDataPointers_d.data(), mDataPointers_d.size());
    ViewBase<std::uint32_t> dataIndices2D[rank-1];
    dataIndices2D[0] = ViewBase<std::uint32_t>(mDataIndices_d.data(), mDataIndices_d.size());

    // host mem
    mModsRedFlat.resize(inFlatSize/K);

    // alloc device mem
    mModsRedFlat_d = std::move(QuICC::Memory::MemBlock<double>(mModsRedFlat.size(), mem.get()));

    // set views
    mModsRedView = modsRed_t(mModsRedFlat_d.data(), mModsRedFlat.size(), modsDims.data(), dataPointers2D, dataIndices2D);
")

set(RadialGpuSetPowerView
"
    mPhysAbsFlat.resize(M * nCols);

    // alloc device mem
    mPhysAbsFlat_d = std::move(QuICC::Memory::MemBlock<double>(mPhysAbsFlat.size(), mem.get()));

    // set device views
    mPhysAbsView = physAbs_t(mPhysAbsFlat_d.data(), mPhysAbsFlat_d.size(), physDims.data(), dataPointers_d, dataIndices_d);
")

set(EnergyCompTree
"
    /// pointwise op on grid points
    mPnt->apply(mPhysAbsView, mPhysView);
    /// Apply reductor to modal space
    mRed->apply(mModsRedView, mPhysAbsView);
")

set(EnergyViewToEigen
"
    assert(static_cast<std::size_t>(rOut.size()) == mModsRedFlat.size());

    for(std::uint32_t i = 0; i < rOut.size(); i++)
    {
        rOut.data()[i]= mModsRedFlat[i];
    }
"
)

set(GpuHostToDevice
"
    cudaErrChk(cudaMemcpy(mModsFlat_d.data(), mModsFlat.data(),
        mModsFlat.size()*sizeof(cuda::std::complex<double>), cudaMemcpyHostToDevice));
"
)

set(GpuPowerDeviceToHost
"
    /// copy back
    cudaErrChk(cudaMemcpy(mModsAbsFlat.data(), mModsAbsFlat_d.data(),
        mModsAbsFlat_d.size()*sizeof(double), cudaMemcpyDeviceToHost));
"
)

set(GpuEnergyDeviceToHost
"
    /// copy back
    cudaErrChk(cudaMemcpy(mModsRedFlat.data(), mModsRedFlat_d.data(),
        mModsRedFlat_d.size()*sizeof(double), cudaMemcpyDeviceToHost));
"
)

set(GpuRadialDeviceToHost
"
    /// copy back
    cudaErrChk(cudaMemcpy(mPhysAbsFlat.data(), mPhysAbsFlat_d.data(),
        mPhysAbsFlat_d.size()*sizeof(double), cudaMemcpyDeviceToHost));
"
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
        set(Grid "${${Kind}Grid}")
        # set(OpInOut "${${Kind}OpInOut}")
        set(SizeRow "${${Kind}SizeRow}")
        set(SizeCol "${${Kind}SizeCol}")
        set(ViewToEigen "${${Kind}ViewToEigen}")
        set(HostToDevice "${${Backend}HostToDevice}")
        set(DeviceToHost "${${Backend}${Kind}DeviceToHost}")

        set(SetReductorView "${${Kind}${Backend}SetReductorView}")
        set(CompTree "${${Kind}CompTree}")

        foreach(Op IN LISTS Ops)

            set(SetIntOp "${IntOp}")
            # manual sub
            string(REGEX REPLACE "@Op@" "${Op}" SetIntOp "${SetIntOp}")
            string(REGEX REPLACE "@BackendNamespace@" "${BackendNamespace}" SetIntOp "${SetIntOp}")

            configure_file(
                "OpV.hpp.in"
                "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/FiniteDiff/Sphere/Reductor/${Op}${Tag}.hpp.inc"
            )
            set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
                PROPERTIES
                PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/FiniteDiff/Sphere/Reductor/${Op}${Tag}.hpp.inc"
            )

            configure_file(
                "OpV.cpp.in"
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/FiniteDiff/Sphere/Reductor/${Op}${Tag}.cpp"
            )
            target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
              "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/FiniteDiff/Sphere/Reductor/${Op}${Tag}.cpp"
            )

        endforeach()
    endforeach()
endforeach()
