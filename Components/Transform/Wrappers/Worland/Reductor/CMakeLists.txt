set(Backends "Cpu")

if(TARGET QuICC::Cuda)
    list(APPEND Backends "Gpu")
endif()

set(Kinds "Power;Energy;Radial;Spectrum")

set(PowerOps "Power;PowerR2;PowerD1R1;PowerSLaplR2")
set(EnergyOps "Energy;EnergyR2;EnergyD1R1;EnergySLaplR2")
set(RadialOps "RadialPower;RadialPowerDivR1;RadialPowerDivR1D1R1")
set(SpectrumOps "Spectrum")

set(PowerOpAttribute "CSL3D")
set(EnergyOpAttribute "CSL3D")
set(RadialOpAttribute "CSL3D")
set(SpectrumOpAttribute "CSL3D")

set(PowerRedType "mods")
set(EnergyRedType "mods")
set(RadialRedType "phys")
set(SpectrumRedType "mods")

set(PowerSizeRow "this->mspSetup->fastSize(0);")
set(PowerSizeCol "this->mspSetup->blockSize();")
set(EnergySizeRow "this->mspSetup->blockSize();")
set(EnergySizeCol "1;")
set(RadialSizeRow "this->mspSetup->fwdSize();")
set(RadialSizeCol "this->mspSetup->blockSize();")
set(SpectrumSizeRow "this->mspSetup->fastSize(0);")
set(SpectrumSizeCol "this->mspSetup->blockSize();")

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

    mInt = std::make_unique<derivedInt_t>(mem);
    derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mInt);
    derivedOp.allocOp(dimensionsInt, layers);
    auto intView = derivedOp.getOp();

    OpsBuilder<op_t, @Op@_t, fwd_t> tBuilderFwd;
    tBuilderFwd.compute(intView, igrid, iweights);
")

set(IntOpAnelastic
"
    // dim 0 - R  - radial modes
    // dim 1 - Nr - radial points
    // dim 2 - L  - harmonic degree
    std::array<std::uint32_t, rank> dimensionsInt {Kanel, Manel, P};

    // Make view integrator operator
    using backendInt_t = @BackendNamespace@::ImplOp<mods_t, phys_t, op_t>;
    using derivedInt_t = Op<mods_t, phys_t, op_t, backendInt_t>;

    mIntAnelastic = std::make_unique<derivedInt_t>(mem);
    derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mIntAnelastic);
    derivedOp.allocOp(dimensionsInt, layers);
    auto intView = derivedOp.getOp();

    OpsBuilder<op_t, @Op@_t, fwd_t> tBuilderFwd;
    tBuilderFwd.compute(intView, igridAnelastic, iweightsAnelastic);
")

# cpu/gpu
set(CpuMemRes "Cpu::NewDelete")
set(GpuMemRes "Cuda::Malloc")

set(GpuLayout "#define QUICC_JW_ALL_ROW_MAJOR")

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

set(CpuSetDataViewsAnelastic
"// set pointers
    ViewBase<std::uint32_t> dataPointers[rank];
    dataPointers[1] = ViewBase<std::uint32_t>(mDataPointersAnelastic.data(), mDataPointersAnelastic.size());
    ViewBase<std::uint32_t> dataIndices[rank];
    dataIndices[1] = ViewBase<std::uint32_t>(mDataIndicesAnelastic.data(), mDataIndicesAnelastic.size());

    // set views
    mModsViewAnelastic = mods_t(mModsFlatAnelastic.data(), mModsFlatAnelastic.size(), modsDims.data(), dataPointers, dataIndices);
    mPhysViewAnelastic = phys_t(mPhysFlatAnelastic.data(), mPhysFlatAnelastic.size(), physDims.data(), dataPointers, dataIndices);"
)

set(PowerCpuSetPowerView
"
    mModsAbsFlat.resize(inFlatSize);
    mModsAbsView = modsAbs_t(mModsAbsFlat.data(), mModsAbsFlat.size(), modsDims.data(), dataPointers, dataIndices);
")

set(PowerCpuSetPowerViewAnelastic
"
    mModsAbsFlatAnelastic.resize(inFlatSize);
    mModsAbsViewAnelastic = modsAbs_t(mModsAbsFlatAnelastic.data(), mModsAbsFlatAnelastic.size(), modsDims.data(), dataPointers, dataIndices);
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

set(EnergyCpuSetReductorViewAnelastic
"
    mModsAbsFlatAnelastic.resize(inFlatSize);
    mModsAbsViewAnelastic = modsAbs_t(mModsAbsFlatAnelastic.data(), mModsAbsFlatAnelastic.size(), modsDims.data(), dataPointers, dataIndices);

    // setup reductor
    mRedAnelastic = std::make_unique<QuICC::Reduction::Cpu::Op<modsRed_t, modsAbs_t, 0>>();

    // set pointers
    ViewBase<std::uint32_t> dataPointers2D[rank-1];
    dataPointers2D[0] = ViewBase<std::uint32_t>(mDataPointersAnelastic.data(), mDataPointersAnelastic.size());
    ViewBase<std::uint32_t> dataIndices2D[rank-1];
    dataIndices2D[0] = ViewBase<std::uint32_t>(mDataIndicesAnelastic.data(), mDataIndicesAnelastic.size());

    // set views
    mModsRedFlatAnelastic.resize(inFlatSize/Kanel);
    mModsRedViewAnelastic = modsRed_t(mModsRedFlatAnelastic.data(), mModsRedFlatAnelastic.size(), modsDims.data(), dataPointers2D, dataIndices2D);"
)

set(RadialCpuSetPowerView
"
    mPhysAbsFlat.resize(M * nCols);
    mPhysAbsView = physAbs_t(mPhysAbsFlat.data(), mPhysAbsFlat.size(), physDims.data(), dataPointers, dataIndices);
")

set(SpectrumCpuSetPowerView
"
    mModsAbsFlat.resize(inFlatSize);
    mModsAbsView = modsAbs_t(mModsAbsFlat.data(), mModsAbsFlat.size(), modsDims.data(), dataPointers, dataIndices);
")

set(SpectrumCpuSetPowerViewAnelastic
"
    mModsAbsFlatAnelastic.resize(inFlatSize);
    mModsAbsViewAnelastic = modsAbs_t(mModsAbsFlatAnelastic.data(), mModsAbsFlatAnelastic.size(), modsDims.data(), dataPointers, dataIndices);
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

set(SpectrumGpuSetPowerView
"
    mModsAbsFlat.resize(inFlatSize);

    // alloc device mem
    mModsAbsFlat_d = std::move(QuICC::Memory::MemBlock<double>(mModsAbsFlat.size(), mem.get()));

    // set device views
    mModsAbsView = modsAbs_t(mModsAbsFlat_d.data(), mModsAbsFlat_d.size(), modsDims.data(), dataPointers_d, dataIndices_d);
")

set(PowerCompTree
"
    if (pF)
    {
        throw std::logic_error(\"Anelastic version of power calculation not fully set up.\");
    }
    /// project to physical space
    mPrj->apply(mPhysView, mModsView);
    /// project to modal space
    mInt->apply(mModsView, mPhysView);
    /// pointwise op on coefficients
    mPnt->apply(mModsAbsView, mModsView);
")

set(EnergyCompTree
"
    /// project to physical space
    if (pF)
    {
        mPrjAnelastic->apply(mPhysViewAnelastic, mModsViewAnelastic);
    }
    else
    {
        mPrj->apply(mPhysView, mModsView);
    }

    if (pF) 
    {
        if (mDensityScale.empty())
        {
            throw std::logic_error(\"mDensityScale is Empty\");
        }
        else
        {
            const std::size_t nR = mPhysViewAnelastic.dims()[0];
            for (std::size_t k = 0; k < mPhysFlatAnelastic.size(); ++k) 
            {
                mPhysFlatAnelastic[k] *= mDensityScale[k % nR];
            }
        }
    }

    if (pF)
    {
        /// project to modal space
        mIntAnelastic->apply(mModsViewAnelastic, mPhysViewAnelastic);
        /// pointwise op on coefficients
        mPntAnelastic->apply(mModsAbsViewAnelastic, mModsViewAnelastic);
        /// reduction op on coefficients
        mRedAnelastic->apply(mModsRedViewAnelastic, mModsAbsViewAnelastic);
    }
    else
    {
        /// project to modal space
        mInt->apply(mModsView, mPhysView);
        /// pointwise op on coefficients
        mPnt->apply(mModsAbsView, mModsView);
        /// reduction op on coefficients
        mRed->apply(mModsRedView, mModsAbsView);
    }
")

set(RadialCompTree
"
    if (pF)
    {
        throw std::logic_error(\"Anelastic version of radial calculation not fully set up.\");
    }
    /// project to physical space
    mPrj->apply(mPhysView, mModsView);
    /// pointwise op on grid points
    mPnt->apply(mPhysAbsView, mPhysView);
")

set(SpectrumCompTree
"
    if (pF) 
    {
        /// project to physical space
        mPrjAnelastic->apply(mPhysViewAnelastic, mModsViewAnelastic);

        if (mDensityScale.empty())
        {
            throw std::logic_error(\"mDensityScale is Empty\");
        }
        else
        {
            const std::size_t nR = mPhysViewAnelastic.dims()[0];
            for (std::size_t k = 0; k < mPhysFlatAnelastic.size(); ++k) 
            {
                mPhysFlatAnelastic[k] *= mDensityScale[k % nR];
            }
        }
    }

    if (pF)
    {
        /// project to modal space
        mIntAnelastic->apply(mModsViewAnelastic, mPhysViewAnelastic);
        /// pointwise op on coefficients
        mPntAnelastic->apply(mModsAbsViewAnelastic, mModsViewAnelastic);
    }
    else
    {
        /// pointwise op on coefficients
        mPnt->apply(mModsAbsView, mModsView);
    }
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
            #ifdef QUICC_JW_ALL_ROW_MAJOR
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
    if (pF)
    {
        assert(static_cast<std::size_t>(rOut.size()) == mModsRedFlatAnelastic.size());

        for (std::uint32_t i = 0; i < rOut.size(); ++i)
        {
            rOut.data()[i] = mModsRedFlatAnelastic[i];
        }
    }
    else
    {
        assert(static_cast<std::size_t>(rOut.size()) == mModsRedFlat.size());

        for (std::uint32_t i = 0; i < rOut.size(); ++i)
        {
            rOut.data()[i] = mModsRedFlat[i];
        }
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
            #ifdef QUICC_JW_ALL_ROW_MAJOR
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

set(SpectrumViewToEigen
"
    if(pF)
    {
        assert(static_cast<std::size_t>(rOut.rows()) <= mModsAbsViewAnelastic.dims()[0]);
        assert(static_cast<std::size_t>(rOut.cols()) <= mModsAbsViewAnelastic.indices()[1].size());
    }
    else
    {
        assert(static_cast<std::size_t>(rOut.rows()) == mModsAbsView.dims()[0]);
        assert(static_cast<std::size_t>(rOut.cols()) == mModsAbsView.indices()[1].size());
    }

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
            #ifdef QUICC_JW_ALL_ROW_MAJOR
                // copy to padded from flattened column and transpose
                if (pF) 
                {
                    outB.data()[i+j*rOut.rows()]= mModsAbsFlatAnelastic[offSet + i*outB.cols()+j];
                }
                else
                {
                    outB.data()[i+j*rOut.rows()]= mModsAbsFlat[offSet + i*outB.cols()+j];
                }
                #else
                // copy to padded from flattened column
                if (pF)
                {
                    outB.data()[i+j*rOut.rows()]= mModsAbsFlatAnelastic[offSet + i+j*(outB.rows()+nf)];
                }
                else
                {
                    outB.data()[i+j*rOut.rows()]= mModsAbsFlat[offSet + i+j*outB.rows()];
                }
                #endif
            }
        }

        if (pF)
        {
            offSet += static_cast<std::int64_t>(outRows + nf) * cols;
        }
        else
        {
            offSet += outB.size();
        }
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

set(GpuSpectrumDeviceToHost
"
    /// copy back
    cudaErrChk(cudaMemcpy(mModsAbsFlat.data(), mModsAbsFlat_d.data(),
        mModsAbsFlat_d.size()*sizeof(double), cudaMemcpyDeviceToHost));
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
    set(SetDataViewsAnelastic "${${Backend}SetDataViewsAnelastic}")
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
        set(SetPowerViewAnelastic "${${Kind}${Backend}SetPowerViewAnelastic}")
        set(SetReductorView "${${Kind}${Backend}SetReductorView}")
        set(SetReductorViewAnelastic "${${Kind}${Backend}SetReductorViewAnelastic}")
        set(SetRadialView "${${Kind}${Backend}SetRadialView}")
        set(SetRadialViewAnelastic "${${Kind}${Backend}SetRadialViewAnelastic}")
        set(CompTree "${${Kind}CompTree}")

        foreach(Op IN LISTS Ops)

            set(SetIntOp "${IntOp}")
            # manual sub
            string(REGEX REPLACE "@Op@" "${Op}" SetIntOp "${SetIntOp}")
            string(REGEX REPLACE "@BackendNamespace@" "${BackendNamespace}" SetIntOp "${SetIntOp}")
            
            set(SetIntOpAnelastic "${IntOpAnelastic}")
            # manual sub
            string(REGEX REPLACE "@Op@" "${Op}" SetIntOpAnelastic "${SetIntOpAnelastic}")
            string(REGEX REPLACE "@BackendNamespace@" "${BackendNamespace}" SetIntOpAnelastic "${SetIntOpAnelastic}")
            
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
