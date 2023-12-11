# Generate wrappers for View operators
set(Kinds "Projector;Integrator")
set(Spaces "Mixed;Complex")
# set(ProjectorMixedOps "P")
# set(ProjectorComplexOps "P")
set(ProjectorMixedOps "P;D1;D2;D3")
set(IntegratorMixedOps "P;D1;D1_P;D1_Neg;D2")
set(ProjectorComplexOps "P;D1;D2;D3;D4;Lapl2D;Df1Lapl2D;Ds1Lapl2D;Mean")
set(IntegratorComplexOps "P;P_Clean;D1;D1_P;D1_Neg;D2;Lapl2D;InvLapl2D;Df1InvLapl2D;Mean")

set(Tags "viewCpu_t")
if(QUICC_USE_CUFFT)
    list(APPEND Tags "viewGpu_t")
endif()
if(QUICC_USE_VKFFT)
    list(APPEND Tags "viewGpuVkFFT_t")
endif()

# Memres
set(viewCpu_tMemRes "Cpu::NewDelete")
set(viewGpu_tMemRes "Cuda::Malloc")
set(viewGpuVkFFT_tMemRes "Cuda::Malloc")
set(ProjectorMem "mem")

# Old interface output type
set(ProjectorMixedOutMatrixType "Matrix")
set(ProjectorMixedInMatrixType "MatrixZ")
set(IntegratorMixedOutMatrixType "MatrixZ")
set(IntegratorMixedInMatrixType "Matrix")
set(ProjectorComplexInMatrixType "MatrixZ")
set(ProjectorComplexOutMatrixType "MatrixZ")
set(IntegratorComplexInMatrixType "MatrixZ")
set(IntegratorComplexOutMatrixType "MatrixZ")

set(ComplexPhysMatrixType "MatrixZ")
set(MixedPhysMatrixType "Matrix")

# New interface types
set(ProjectorArgInType "mods_t")
set(ProjectorArgOutType "phys_t")
set(IntegratorArgInType "phys_t")
set(IntegratorArgOutType "mods_t")

# Copy to padded
set(viewCpu_tProjectorMixedCopyPad "// dealias or overintegration
    std::size_t nDealias = mspSetup->specSize();
    mPadded.topRows(nDealias) = in.topRows(nDealias);")
set(viewGpu_tProjectorMixedCopyPad "// dealias or overintegration
    std::size_t nDealias = mspSetup->specSize();

    // cpu -> gpu data
    auto dSz = sizeof(decltype(mMods_d)::type);
    cudaErrChk(cudaMemcpy2DAsync(mMods_d.data(), mPadded.rows()*dSz, in.data(), in.rows()*dSz,
        nDealias*dSz, in.cols(), cudaMemcpyHostToDevice));")
set(viewGpuVkFFT_tProjectorMixedCopyPad "${viewGpu_tProjectorMixedCopyPad}")
set(viewGpu_tIntegratorMixedCopyPad "// cpu -> gpu data
    auto dSz = sizeof(decltype(mPhys_d)::type);
    cudaErrChk(cudaMemcpyAsync(mPhys_d.data(), in.data(), mPhys_d.size()*dSz,
        cudaMemcpyHostToDevice));")
set(viewGpuVkFFT_tIntegratorMixedCopyPad "${viewGpu_tIntegratorMixedCopyPad}")
set(viewCpu_tProjectorComplexCopyPad "// dealias or overintegration
    std::size_t nDealias = mspSetup->specSize();
    auto negN = nDealias / 2;
    auto posN = negN + nDealias % 2;
    mPadded.topRows(posN) = in.topRows(posN);
    mPadded.bottomRows(negN) = in.bottomRows(negN);")
set(viewGpu_tProjectorComplexCopyPad "${viewCpu_tProjectorComplexCopyPad}

    // cpu -> gpu data
    auto dSz = sizeof(decltype(mMods_d)::type);
    cudaErrChk(cudaMemcpyAsync(mMods_d.data(), mPadded.data(), mMods_d.size()*dSz,
        cudaMemcpyHostToDevice));"
    )
set(viewGpuVkFFT_tProjectorComplexCopyPad "${viewGpu_tProjectorComplexCopyPad}")
set(viewGpu_tIntegratorComplexCopyPad "// cpu -> gpu data
    auto dSz = sizeof(decltype(mPhys_d)::type);
    cudaErrChk(cudaMemcpyAsync(mPhys_d.data(), in.data(), mPhys_d.size()*dSz,
        cudaMemcpyHostToDevice));"
    )
set(viewGpuVkFFT_tIntegratorComplexCopyPad "${viewGpu_tIntegratorComplexCopyPad}")

set(MixedViewAttribute "DCCSC3D")
set(ComplexViewAttribute "DCCSC3DInOrder")

# P - Projector
set(ProjectorPviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
Cpu::DiffOp<mods_t, mods_t, 0, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorPviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorPviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, bwd_t>>")
endif()
# P - Integrator
set(IntegratorPviewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 0, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorPviewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorPviewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, fwd_t>>")
endif()
# P_Clean - Integrator
set(IntegratorP_CleanviewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 0, fwd_t, zeroResetMean_m>>")
if(QUICC_USE_VKFFT)
set(IntegratorP_CleanviewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, fwd_t, zeroResetMean_m>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorP_CleanviewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 0, fwd_t, zeroResetMean_m>>")
endif()
# D1 - Projector
set(ProjectorD1viewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
Cpu::DiffOp<mods_t, mods_t, 1, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorD1viewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorD1viewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, bwd_t>>")
endif()
# D1 - Integrator
set(IntegratorD1viewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 1, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorD1viewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorD1viewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>")
endif()
# D1_P - Integrator
set(IntegratorD1_PviewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_m>>")
if(QUICC_USE_VKFFT)
set(IntegratorD1_PviewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_m>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorD1_PviewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_m>>")
endif()
# D1_Neg - Integrator
set(IntegratorD1_NegviewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_m>>")
if(QUICC_USE_VKFFT)
set(IntegratorD1_NegviewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_m>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorD1_NegviewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_m>>")
endif()
# D2 - Projector
set(ProjectorD2viewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::DiffOp<mods_t, mods_t, 2, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorD2viewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 2, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorD2viewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 2, bwd_t>>")
endif()
# D2 - Integrator
set(IntegratorD2viewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 2, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorD2viewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 2, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorD2viewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 2, fwd_t>>")
endif()
# D3 - Projector
set(ProjectorD3viewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::DiffOp<mods_t, mods_t, 3, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorD3viewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 3, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorD3viewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 3, bwd_t>>")
endif()
# D3 - Integrator
set(IntegratorD3viewCpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
            Cpu::DiffOp<mods_t, mods_t, 3, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorD3viewGpuVkFFT_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::VkFft::FftOp<<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 3, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorD3viewGpu_tPtrParam "DOp<mods_t, phys_t,
        QuICC::Fft::CuFft::FftOp<<mods_t, phys_t>,
            Cuda::DiffOp<mods_t, mods_t, 3, fwd_t>>")
endif()
# D4 - Projector
set(ProjectorD4viewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::DiffOp<mods_t, mods_t, 4, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorD4viewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 4, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorD4viewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::DiffOp<mods_t, mods_t, 4, bwd_t>>")
endif()
# Lapl2D - Projector
set(ProjectorLapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorLapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorLapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t>>")
endif()
# Lapl2D - Integrator
set(IntegratorLapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorLapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorLapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>>")
endif()
# InvLapl2D - Integrator
set(IntegratorInvLapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>>")
if(QUICC_USE_VKFFT)
set(IntegratorInvLapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorInvLapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>>")
endif()
# Df1InvLapl2D - Integrator
set(IntegratorDf1InvLapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>,
            Cpu::DiffOp<mods_t, mods_t, 1, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorDf1InvLapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorDf1InvLapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>,
            Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>")
endif()
# Df1Lapl2D - Projector
set(ProjectorDf1Lapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 3, 0, 1, 2, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorDf1Lapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 3, 0, 1, 2, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorDf1Lapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 3, 0, 1, 2, bwd_t>>")
endif()
# Ds1Lapl2D - Projector
set(ProjectorDs1Lapl2DviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::Diff2DOp<mods_t, mods_t, 2, 1, 0, 3, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorDs1Lapl2DviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 1, 0, 3, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorDs1Lapl2DviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::Diff2DOp<mods_t, mods_t, 2, 1, 0, 3, bwd_t>>")
endif()
# Mean - Projector
set(ProjectorMeanviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::MeanOp<mods_t, mods_t, bwd_t>>")
if(QUICC_USE_VKFFT)
set(ProjectorMeanviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::MeanOp<mods_t, mods_t, bwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(ProjectorMeanviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::MeanOp<mods_t, mods_t, bwd_t>>")
endif()
# Mean - Integrator
set(IntegratorMeanviewCpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
            Cpu::MeanOp<mods_t, mods_t, fwd_t>>")
if(QUICC_USE_VKFFT)
set(IntegratorMeanviewGpuVkFFT_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::VkFft::FftOp<phys_t, mods_t>,
            Cuda::MeanOp<mods_t, mods_t, fwd_t>>")
endif()
if(QUICC_USE_CUFFT)
set(IntegratorMeanviewGpu_tPtrParam "DOp<phys_t, mods_t,
        QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
            Cuda::MeanOp<mods_t, mods_t, fwd_t>>")
endif()

# View Cpu - Projector
set(ProjectorviewCpu_tViewStorage "// views
    mods_t mods({const_cast<std::remove_reference_t<decltype(in)>::Scalar*>(mPadded.data()),
        static_cast<std::size_t>(mPadded.size())},
        {{static_cast<std::uint32_t>(nDealias),
            columns, layers}}, pointers, indices, mPadded.rows());
    phys_t phys({rOut.data(),
        static_cast<std::size_t>(rOut.size())},
            {{static_cast<std::uint32_t>(rOut.rows()),
            columns, layers}}, pointers, indices);")

# View Cpu - Integrator
set(IntegratorviewCpu_tViewStorage "// views
    mods_t mods({rOut.data(),
        static_cast<std::size_t>(rOut.size())},
        {{static_cast<std::uint32_t>(rOut.rows()),
            columns, layers}}, pointers, indices);
    phys_t phys({const_cast<std::remove_reference_t<decltype(in)>::Scalar*>(in.data()),
        static_cast<std::size_t>(in.size())},
            {{static_cast<std::uint32_t>(in.rows()),
            columns, layers}}, pointers, indices);")


set(ProjectorviewGpu_tViewStorage "
    // device mem for metadata
    QuICC::Memory::MemBlock<std::uint32_t> ptrAndIdx_d(psz + isz, mMemRes.get());

    // Set device pointers
    constexpr size_t rank = 3;
    ViewBase<std::uint32_t> pointers_d[rank];
    pointers_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data(), psz);
    ViewBase<std::uint32_t> indices_d[rank];
    indices_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data() + psz, isz);

    // dims
    std::array<std::uint32_t, rank> dimMods {{static_cast<std::uint32_t>(mPadded.rows()),
        columns, layers}};
    std::array<std::uint32_t, rank> dimPhys {{static_cast<std::uint32_t>(rOut.rows()),
        columns, layers}};

    // views
    mods_t mods(reinterpret_cast<MatrixZ::Scalar*>(mMods_d.data()), mMods_d.size(),
        dimMods.data(), pointers_d, indices_d);
    phys_t phys(reinterpret_cast<std::remove_reference_t<decltype(rOut)>::Scalar*>(mPhys_d.data()), mPhys_d.size(),
        dimPhys.data(), pointers_d, indices_d);

    // cpu -> gpu metadata
    cudaErrChk(cudaMemcpyAsync(mods.pointers()[1].data(), pointers[1].data(),
        psz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpyAsync(mods.indices()[1].data(), indices[1].data(),
        isz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));"
)
set(ProjectorviewGpuVkFFT_tViewStorage ${ProjectorviewGpu_tViewStorage})

set(IntegratorviewGpu_tViewStorage "
    // device mem for metadata
    QuICC::Memory::MemBlock<std::uint32_t> ptrAndIdx_d(psz + isz, mMemRes.get());

    // Set device pointers
    constexpr size_t rank = 3;
    ViewBase<std::uint32_t> pointers_d[rank];
    pointers_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data(), psz);
    ViewBase<std::uint32_t> indices_d[rank];
    indices_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data() + psz, isz);

    // dims
    std::array<std::uint32_t, rank> dimMods {{static_cast<std::uint32_t>(rOut.rows()),
        columns, layers}};
    std::array<std::uint32_t, rank> dimPhys {{static_cast<std::uint32_t>(in.rows()),
        columns, layers}};

    // views
    mods_t mods(reinterpret_cast<MatrixZ::Scalar*>(mMods_d.data()), mMods_d.size(),
        dimMods.data(), pointers_d, indices_d);
    phys_t phys(reinterpret_cast<std::remove_reference_t<decltype(in)>::Scalar*>(mPhys_d.data()), mPhys_d.size(),
        dimPhys.data(), pointers_d, indices_d);

    // cpu -> gpu metadata
    cudaErrChk(cudaMemcpyAsync(mods.pointers()[1].data(), pointers[1].data(),
        psz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpyAsync(mods.indices()[1].data(), indices[1].data(),
        isz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));"
)
set(IntegratorviewGpuVkFFT_tViewStorage ${IntegratorviewGpu_tViewStorage})

set(ProjectorviewGpu_tCopyBack "// gpu -> cpu
    cudaErrChk(cudaMemcpy(rOut.data(), mPhys_d.data(),
        mPhys_d.size()*sizeof(decltype(mPhys_d)::type), cudaMemcpyDeviceToHost));")
set(ProjectorviewGpuVkFFT_tCopyBack ${ProjectorviewGpu_tCopyBack})


set(IntegratorviewGpu_tCopyBack "// gpu -> cpu
    cudaErrChk(cudaMemcpy(rOut.data(), mMods_d.data(),
        mMods_d.size()*sizeof(decltype(mMods_d)::type), cudaMemcpyDeviceToHost));")
set(IntegratorviewGpuVkFFT_tCopyBack ${IntegratorviewGpu_tCopyBack})

# Indices/Pointers
set(ComplexIndicesPointers "// get 3D dimensions
    auto idBlocks = mspSetup->idBlocks();
    std::uint32_t activeLayers = idBlocks.rows();

    // here we guestimate based on highest index
    // since we don't how how the problem is setup globally
    std::uint32_t layers = idBlocks(activeLayers-1,0)+1;

    std::uint32_t columns = 0;
    for(std::uint32_t i = 0; i < activeLayers; ++i)
    {
        std::uint32_t maxCol = idBlocks(i, 1);
        columns = std::max(maxCol, columns);
    }

    std::uint32_t activeColumns = mPadded.cols();

    std::array<std::vector<std::uint32_t>, 3> pointers{};
    std::array<std::vector<std::uint32_t>, 3> indices{};
    auto psz = layers+1;
    auto isz = activeColumns;
    pointers[1].resize(psz);
    indices[1].resize(isz);

    // set pointers and indices
    pointers[1][0] = 0;
    std::uint32_t idk = 0;
    for(std::uint32_t ii = 0; ii < layers; ++ii)
    {
        std::uint32_t k = idBlocks(idk, 0);
        std::uint32_t numCol = 0;
        if (k == ii)
        {
            numCol = idBlocks(idk, 1);
            idk += 1;
        }

        pointers[1][ii+1] = pointers[1][ii] + numCol;

        std::uint32_t j = 0;
        for(std::uint32_t idj = pointers[1][ii]; idj < pointers[1][ii+1]; ++idj)
        {
            // making j index up
            indices[1][idj] = j;
            j += 1;
        }
    }"
)

set(MixedIndicesPointers "// get 3D dimensions
    // here we make things up
    // since we don't how how the problem is setup globally

    std::uint32_t layers = 1;
    std::uint32_t columns = mPadded.cols();

    std::array<std::vector<std::uint32_t>, 3> pointers{};
    std::array<std::vector<std::uint32_t>, 3> indices{};
    auto psz = layers+1;
    auto isz = columns;
    pointers[1].resize(psz);
    indices[1].resize(isz);

    // set pointers and indices
    pointers[1][0] = 0;
    pointers[1][1] = columns;
    for(std::uint32_t j = 0; j < columns; ++j)
    {
        // making j index up
        indices[1][j] = j;
    }"
)

# PreOp
set(ProjectorPreOp "/**
* @brief Apply pre FFT operator
*
* @param in   Input mods values
* @param out  Copied and padded input
*/
void applyPreOperator(MatrixZ& out, const MatrixZ& in) const final {};"
)

# PostOp
set(IntegratorPostOp "/**
* @brief Apply post FFT operator
*
* @param rOut Output values
*/
void applyPostOperator(MatrixZ& rOut) const final {};"
)

# loop
foreach(Space IN LISTS Spaces)
    set(ViewAttribute "${${Space}ViewAttribute}")
    foreach(Kind IN LISTS Kinds)
        set(Ops "${${Kind}${Space}Ops}")
        set(Mem "${${Kind}Mem}")
        set(OutMatrixType "${${Kind}${Space}OutMatrixType}")
        set(InMatrixType "${${Kind}${Space}InMatrixType}")

        set(ArgInType "${${Kind}ArgInType}")
        set(ArgOutType "${${Kind}ArgOutType}")
        set(ArgsType "${ArgOutType}, ${ArgInType}")
        string(REGEX REPLACE "_t" "" Args "${ArgsType}")

        set(PreOp "${${Kind}PreOp}")
        set(PostOp "${${Kind}PostOp}")

        # Complex headers
        set(ComplexCppHeader "#include \"ViewOps/Fourier/${Space}/Diff2D.hpp\"
#include \"ViewOps/Fourier/${Space}/Mean.hpp\""
        )

        # Cuda only bits
        set(viewGpu_tCuhHeader "#include <cuComplex.h>
#include \"Memory/Memory.hpp\"
#include \"Memory/Cuda/Malloc.hpp\""
        )
        set(viewGpuVkFFT_tCuhHeader "${viewGpu_tCuhHeader}")

        set(viewGpu_tMemBlock "/// Device memory
    mutable std::unique_ptr<QuICC::Memory::memory_resource> mMemRes;
    mutable QuICC::Memory::MemBlock<cuDoubleComplex> mMods_d;
    mutable QuICC::Memory::MemBlock<${${Space}PhysMatrixType}::Scalar> mPhys_d;")
        set(viewGpuVkFFT_tMemBlock "${viewGpu_tMemBlock}")

        set(viewGpu_tCuHeader "#include <cuda_runtime_api.h>
    #include \"Memory/Cuda/Malloc.hpp\"
    #include \"Cuda/CudaUtil.hpp\""
        )
        set(viewGpuVkFFT_tCuHeader "${viewGpu_tCuHeader}")

        set(viewGpu_tAlloc "// Alloc device memory
    mMemRes = std::make_unique<QuICC::Memory::Cuda::Malloc>();
    mMods_d = std::move(QuICC::Memory::MemBlock<cuDoubleComplex>(bwdSize*blockSize, mMemRes.get()));
    auto fwdSize = mspSetup->fwdSize();
    mPhys_d = std::move(QuICC::Memory::MemBlock<${${Space}PhysMatrixType}::Scalar>(fwdSize*blockSize, mMemRes.get()));

    // zero gpu mem
    cudaErrChk(cudaMemcpy(mMods_d.data(), mPadded.data(),
        mMods_d.size()*sizeof(decltype(mMods_d)::type), cudaMemcpyHostToDevice));
        "
        )
        set(viewGpuVkFFT_tAlloc "${viewGpu_tAlloc}")

            foreach(Op IN LISTS Ops)
            foreach(Tag IN LISTS Tags)
                set(CopyPad "${${Tag}${Kind}${Space}CopyPad}")
                set(CuhHeader "${${Tag}CuhHeader}")
                set(MemBlock "${${Tag}MemBlock}")
                set(MemRes "${${Tag}MemRes}")
                set(IndicesPointers "${${Space}IndicesPointers}")
                configure_file(
                    "OpV.hpp.in"
                    "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Fourier/${Space}/${Kind}/${Op}${Tag}.hpp.inc"
                )
                set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
                    PROPERTIES
                        PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Fourier/${Space}/${Kind}/${Op}${Tag}.hpp.inc"
                )
                set(CppHeader "${${Space}CppHeader}")
                set(CuHeader "${${Tag}CuHeader}")
                set(PtrParam "${${Kind}${Op}${Tag}PtrParam}")
                set(Alloc "${${Tag}Alloc}")
                set(ViewStorage "${${Kind}${Tag}ViewStorage}")
                set(CopyBack "${${Kind}${Tag}CopyBack}")

                message(DEBUG "CppHeader: ${CppHeader}")
                message(DEBUG "PtrParam: ${PtrParam}")
                configure_file(
                    "OpV.cpp.in"
                    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Fourier/${Space}/${Kind}/${Op}${Tag}.cpp"
                )
                target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
                    "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Fourier/${Space}/${Kind}/${Op}${Tag}.cpp"
                )
            endforeach()
        endforeach()
    endforeach()
endforeach()
