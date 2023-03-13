
#include <iostream>
#include <complex>

#include "D.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Mixed/Diff.hpp"
#include "Fft/Fft.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
namespace Projector {

using namespace QuICC::Memory;

using mods_t = View<std::complex<double>, DCCSC3D>;
using phys_t = View<double, DCCSC3D>;

template<class Tout, class Tin, class FftBackend, class DiffBackend>
DOp<Tout, Tin, FftBackend, DiffBackend>::DOp(ScaleType scale) : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>(scale))
{
}

template<class Tout, class Tin, class FftBackend, class DiffBackend>
DOp<Tout, Tin, FftBackend, DiffBackend>::DOp() : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>())
{
}

template<class Tout, class Tin, class FftBackend, class DiffBackend>
void DOp<Tout, Tin, FftBackend, DiffBackend>::applyImpl(Tout& out, Tin& in)
{
    Profiler::RegionFixture<4> fix("DOp::applyImpl");

    // differentiate in place
    mDiff->apply(in, in);

    // FFT
    mFft->apply(out, in);
}

template<class Tout, class Tin, class FftBackend, class DiffBackend>
void DOp<Tout, Tin, FftBackend, DiffBackend>::applyImpl(Tout& out, const Tin& in)
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
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::DiffOp<mods_t, mods_t, 0, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::DiffOp<mods_t, mods_t, 1, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::DiffOp<mods_t, mods_t, 2, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::DiffOp<mods_t, mods_t, 3, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
    Cpu::DiffOp<mods_t, mods_t, 4, bwd_t>>;
#ifdef QUICC_USE_CUFFT
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::DiffOp<mods_t, mods_t, 0, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::DiffOp<mods_t, mods_t, 1, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::DiffOp<mods_t, mods_t, 2, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::DiffOp<mods_t, mods_t, 3, bwd_t>>;
template class DOp<phys_t, mods_t,
    QuICC::Fft::CuFft::FftOp<phys_t, mods_t>,
    Cuda::DiffOp<mods_t, mods_t, 4, bwd_t>>;
#endif

} // namespace Projector
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
