set(Backends "Cpu")

if(TARGET QuICC::Cuda)
    list(APPEND Backends "Gpu")
endif()

set(Kinds "Power;Energy;Radial")

set(PowerOps "Power;PowerR2;PowerD1R1;PowerSLaplR2")
set(EnergyOps "Energy;EnergyR2;EnergyD1R1;EnergySLaplR2")
set(RadialOps "RadialPower;RadialPowerDivR1;RadialPowerDivR1D1R1")

set(PowerOpAttribute "CSL3D")
set(EnergyOpAttribute "CSL3D")
set(RadialOpAttribute "CSL3D")

set(PowerRedType "mods")
set(EnergyRedType "mods")
set(RadialRedType "phys")

set(PowerGrid "this->computePowerQuadrature(igrid, iweights, icompgrid.size());")
set(EnergyGrid "this->computePowerQuadrature(igrid, iweights, icompgrid.size());")
set(RadialGrid "igrid = icompgrid; iweights = icompweights;")

set(PowerSizeRow "this->mspSetup->fastSize(0);")
set(PowerSizeCol "this->mspSetup->blockSize();")
set(EnergySizeRow "this->mspSetup->blockSize();")
set(EnergySizeCol "1;")
set(RadialSizeRow "this->mspSetup->fwdSize();")
set(RadialSizeCol "this->mspSetup->blockSize();")

# Ops
set(PowerPolyBuilderBwd "r_1Wnl<QuICC::Polynomial::Worland::recurrence_t>")
set(PowerPolyBuilderFwd "Wnl")
set(PowerMeanBuilderFwd "void")
set(PowerlShift "-1")

set(PowerR2PolyBuilderBwd "Wnl")
set(PowerR2PolyBuilderFwd "Wnl")

set(PowerD1R1PolyBuilderBwd "r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>")
set(PowerD1R1PolyBuilderFwd "Wnl")
set(PowerD1R1MeanBuilderFwd "void")
set(PowerD1R1lShift "-1")

set(PowerSLaplR2PolyBuilderBwd "slaplWnl")
set(PowerSLaplR2PolyBuilderFwd "Wnl")
set(PowerSLaplR2Output "mModsAbsFlat")

set(EnergyPolyBuilderBwd "r_1Wnl<QuICC::Polynomial::Worland::recurrence_t>")
set(EnergyPolyBuilderFwd "Wnl")
set(EnergyMeanBuilderFwd "void")
set(EnergylShift "-1")

set(EnergyR2PolyBuilderBwd "Wnl")
set(EnergyR2PolyBuilderFwd "Wnl")

set(EnergyD1R1PolyBuilderBwd "r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>")
set(EnergyD1R1PolyBuilderFwd "Wnl")
set(EnergyD1R1MeanBuilderFwd "void")
set(EnergyD1R1lShift "-1")

set(EnergySLaplR2PolyBuilderBwd "slaplWnl")
set(EnergySLaplR2PolyBuilderFwd "Wnl")
set(EnergySLaplR2Output "mModsAbsFlat")

set(RadialPowerPolyBuilderBwd "Wnl")

set(RadialPowerDivR1PolyBuilderBwd "r_1Wnl<QuICC::Polynomial::Worland::recurrence_t>")
set(RadialPowerDivR1MeanBuilderBwd "void")

set(RadialPowerDivR1D1R1PolyBuilderBwd "r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>")
set(RadialPowerDivR1D1R1MeanBuilderBwd "void")

# Integrator operator
set(IntOp
"
    // dim 0 - R  - radial modes
    // dim 1 - Nr - radial points
    // dim 2 - L  - harmonic degree
    std::array<std::uint32_t, rank> dimensionsInt {K, M, P};

    // Make view integrator operator
    using backendInt_t = @BackendNamespace@::ImplOp<mods_t, phys_t, op_t>;
    using derivedInt_t = Op<mods_t, phys_t, op_t, backendInt_t>;

    mInt = std::make_unique<derivedInt_t>(dimensionsInt, layers, mem);
    derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mInt);
    auto intView = derivedOp.getOp();

    Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
        Polynomial::Worland::worland_sphenergy_t::DBETA, @lShift@);
    QuICC::DenseSM::Worland::@OpBuilderFwd@ denseBuilder(fWnl);
    Builder<op_t, QuICC::DenseSM::Worland::@OpBuilderFwd@, fwd_t> tBuilderFwd(denseBuilder);
    tBuilderFwd.compute(intView, igrid, iweights);
")

# cpu/gpu
set(CpuMemRes "Cpu::NewDelete")
set(GpuMemRes "Cuda::Malloc")

set(GpuLayout "#define QUICC_JW_ROW_MAJOR")

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
    mPhysView = phys_t(mPhysFlat.data(), mPhysFlat.size(), physDims.data(), dataPointers, dataIndices);"
)

set(PowerCpuSetPowerView
"
    mModsAbsFlat.resize(inFlatSize);
    mModsAbsView = modsAbs_t(mModsAbsFlat.data(), mModsAbsFlat.size(), modsDims.data(), dataPointers, dataIndices);
")

set(EnergyCpuSetReductorView
"
    mModsAbsFlat.resize(inFlatSize);
    mModsAbsView = modsAbs_t(mModsAbsFlat.data(), mModsAbsFlat.size(), modsDims.data(), dataPointers, dataIndices);

    // setup reductor
    mRed = std::make_unique<QuICC::Reduction::Cpu::Op<modsRed_t, modsAbs_t, 0>>();

    // set pointers
    ViewBase<std::uint32_t> dataPointers2D[rank-1];
    dataPointers2D[0] = ViewBase<std::uint32_t>(mDataPointers.data(), mDataPointers.size());
    ViewBase<std::uint32_t> dataIndices2D[rank-1];
    dataIndices2D[0] = ViewBase<std::uint32_t>(mDataIndices.data(), mDataIndices.size());

    // set views
    mModsRedFlat.resize(inFlatSize/K);
    mModsRedView = modsRed_t(mModsRedFlat.data(), mModsRedFlat.size(), modsDims.data(), dataPointers2D, dataIndices2D);"
)

set(RadialCpuSetPowerView
"
    mPhysAbsFlat.resize(M * nCols);
    mPhysAbsView = physAbs_t(mPhysAbsFlat.data(), mPhysAbsFlat.size(), physDims.data(), dataPointers, dataIndices);
")

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

set(PowerCompTree
"
    /// project to modal space
    mInt->apply(mModsView, mPhysView);
    /// pointwise op on coefficients
    mPnt->apply(mModsAbsView, mModsView);
")

set(EnergyCompTree
"
    /// project to modal space
    mInt->apply(mModsView, mPhysView);
    /// pointwise op on coefficients
    mPnt->apply(mModsAbsView, mModsView);
    /// reduction op on coefficients
    mRed->apply(mModsRedView, mModsAbsView);
")

set(RadialCompTree
"
    /// pointwise op on grid points
    mPnt->apply(mPhysAbsView, mPhysView);
")

set(PowerViewToEigen
"
    assert(static_cast<std::size_t>(rOut.rows()) == mModsAbsView.dims()[0]);
    assert(static_cast<std::size_t>(rOut.cols()) == mModsAbsView.indices()[1].size());

    start = 0;
    offSet = 0;
    for(std::uint32_t p = 0; p < nLayers; p++)
    {
        int cols = this->mspSetup->mult(p);
        int outRows = this->mspSetup->fastSize(p);

        Eigen::Ref<Matrix> outB = rOut.block(0, start, outRows, cols);

        for (std::int64_t j = 0; j < outB.cols(); ++j)
        {
            for (std::int64_t i = 0; i < outB.rows(); ++i)
            {
                #ifdef QUICC_JW_ROW_MAJOR
                // copy to padded from flattened column and transpose
                outB.data()[i+j*rOut.rows()]= mModsAbsFlat[offSet + i*outB.cols()+j];
                #else
                // copy to padded from flattened column
                outB.data()[i+j*rOut.rows()]= mModsAbsFlat[offSet + i+j*outB.rows()];
                #endif
            }
        }

        offSet += outB.size();
        start += cols;
    }
"
)

set(EnergyViewToEigen
"
    assert(static_cast<std::size_t>(rOut.size()) == mModsRedFlat.size());

    for(std::uint32_t i = 0; i < rOut.size(); i++)
    {
        rOut.data()[i]= mModsRedFlat[i];
    }
"
)

set(RadialViewToEigen
"
    assert(static_cast<std::size_t>(rOut.rows()) == mPhysAbsView.dims()[0]);
    assert(static_cast<std::size_t>(rOut.cols()) == mPhysAbsView.indices()[1].size());

    start = 0;
    offSet = 0;
    int outRows = this->mspSetup->fwdSize();
    for(std::uint32_t p = 0; p < nLayers; p++)
    {
        int cols = this->mspSetup->mult(p);

        Eigen::Ref<Matrix> outB = rOut.block(0, start, outRows, cols);

        for (std::int64_t j = 0; j < outB.cols(); ++j)
        {
            for (std::int64_t i = 0; i < outB.rows(); ++i)
            {
                #ifdef QUICC_JW_ROW_MAJOR
                // copy to padded from flattened column and transpose
                outB.data()[i+j*rOut.rows()]= mPhysAbsFlat[offSet + i*outB.cols()+j];
                #else
                // copy to padded from flattened column
                outB.data()[i+j*rOut.rows()]= mPhysAbsFlat[offSet + i+j*outB.rows()];
                #endif
            }
        }

        offSet += outB.size();
        start += cols;
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
        set(RedType "${${Kind}RedType}")
        set(Grid "${${Kind}Grid}")
        # set(OpInOut "${${Kind}OpInOut}")
        set(SizeRow "${${Kind}SizeRow}")
        set(SizeCol "${${Kind}SizeCol}")
        set(ViewToEigen "${${Kind}ViewToEigen}")
        set(HostToDevice "${${Backend}HostToDevice}")
        set(DeviceToHost "${${Backend}${Kind}DeviceToHost}")

        set(SetPowerView "${${Kind}${Backend}SetPowerView}")
        set(SetReductorView "${${Kind}${Backend}SetReductorView}")
        set(SetRadialView "${${Kind}${Backend}SetRadialView}")
        set(CompTree "${${Kind}CompTree}")

        foreach(Op IN LISTS Ops)

            set(PolyBuilderFwd "${${Op}PolyBuilderFwd}")
            set(PolyBuilderBwd "${${Op}PolyBuilderBwd}")


            set(OpBuilderFwd "Operator<${PolyBuilderFwd}>")
            set(OpBuilderBwd "Operator<${PolyBuilderBwd}>")
            set(lShift "${${Op}lShift}")
            if(NOT lShift)
                set(lShift "0")
            endif()

            set(MeanBuilderFwd "${${Op}MeanBuilderFwd}")
            if(MeanBuilderFwd)
                if (NOT (MeanBuilderFwd STREQUAL "void"))
                    set(MeanHeader "#include \"QuICC/Polynomial/Worland/${MeanBuilderFwd}.hpp\"")
                endif()
                set(OpBuilderFwd "OperatorWithMean<${PolyBuilderFwd}, ${MeanBuilderFwd}>")
            endif()

            set(MeanBuilderBwd "${${Op}MeanBuilderBwd}")
            if(MeanBuilderBwd)
                if (NOT (MeanBuilderBwd STREQUAL "void"))
                    set(MeanHeader "#include \"QuICC/Polynomial/Worland/${MeanBuilderBwd}.hpp\"")
                endif()
                set(OpBuilderBwd "OperatorWithMean<${PolyBuilderBwd}, ${MeanBuilderBwd}>")
            endif()


            if(PolyBuilderFwd)
                set(SetIntOp "${IntOp}")
                # manual sub
                string(REGEX REPLACE "@OpBuilderFwd@" "${OpBuilderFwd}" SetIntOp "${SetIntOp}")
                string(REGEX REPLACE "@lShift@" "${lShift}" SetIntOp "${SetIntOp}")
                string(REGEX REPLACE "@BackendNamespace@" "${BackendNamespace}" SetIntOp "${SetIntOp}")
            endif()

            set(PolyBuilderScale "${${Op}PolyBuilderScale}")

            configure_file(
                "OpV.hpp.in"
                "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Worland/Reductor/${Op}${Tag}.hpp.inc"
            )
            set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
                PROPERTIES
                    PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Worland/Reductor/${Op}${Tag}.hpp.inc"
            )

            configure_file(
                "OpV.cpp.in"
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Worland/Reductor/${Op}${Tag}.cpp"
            )
            target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
                "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Worland/Reductor/${Op}${Tag}.cpp"
            )

        endforeach()
    endforeach()
endforeach()
