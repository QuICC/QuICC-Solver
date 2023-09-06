#include <iostream>
#include <complex>

#include "D.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Complex/Diff.hpp"
#include "ViewOps/Fourier/Complex/Diff2D.hpp"
#include "ViewOps/Fourier/Complex/Mean.hpp"
#include "Fft/Fft.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
namespace Integrator {

using namespace QuICC::Memory;

using mods_t = View<std::complex<double>, DCCSC3D>;
using phys_t = View<std::complex<double>, DCCSC3D>;

template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::DOp(ScaleType scale) : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>(scale))
{
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2 = std::make_unique<DiffBackend2>(scale);
    }
}

template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::DOp() : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>())
{
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2 = std::make_unique<DiffBackend2>();
    }
}


template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
void DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::applyImpl(Tout& out, Tin& in)
{
    Profiler::RegionFixture<4> fix("DOp::applyImpl");

    // FFT
    mFft->apply(out, in);

    // differentiate in place
    mDiff->apply(out, out);

    // second optional stage needed for Df1InvLapl2D
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2->apply(out, out);
    }
}

template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
void DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::applyImpl(Tout& out, const Tin& in)
{
    // tmp storage
    std::vector<typename Tin::ScalarType> tmp(in.size());

    // copy
    for(std::size_t i = 0; i < in.size(); ++i)
    {
        tmp.data()[i] = in.data()[i];
    }

    // tmp view
    mods_t tmpV(tmp.data(), in.size(), in.dims(), in.pointers(), in.indices());

    this->applyImpl(out, tmpV);
}

// explicit instantations
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 0, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 0, fwd_t, zeroResetMean_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 1, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 2, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 3, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::Fftw::FftOp<mods_t, phys_t>,
    Cpu::DiffOp<mods_t, mods_t, 4, fwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>,
    Cpu::DiffOp<mods_t, mods_t, 1, fwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::MeanOp<mods_t, mods_t, fwd_t>>;

#ifdef QUICC_USE_CUFFT
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 0, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 0, fwd_t, zeroResetMean_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_m>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 2, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 3, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::DiffOp<mods_t, mods_t, 4, fwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>,
    Cuda::DiffOp<mods_t, mods_t, 1, fwd_t>>;
template class DOp<mods_t, phys_t,
    QuICC::Fft::CuFft::FftOp<mods_t, phys_t>,
    Cuda::MeanOp<mods_t, mods_t, fwd_t>>;
#endif

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
