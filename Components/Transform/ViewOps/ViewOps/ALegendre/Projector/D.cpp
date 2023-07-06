
#include <iostream>
#include <complex>

#include "D.hpp"
#include "View/View.hpp"
#include "ViewOps/ALegendre/Diff.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {
namespace Projector {

using namespace QuICC::Memory;

// using mods_t = View<std::complex<double>, DCCSC3D>;
// using phys_t = View<double, DCCSC3D>;

template<class Tout, class Tin, class Backend>
DOp<Tout, Tin, Backend>::DOp(ScaleType scale) :
    mOp(std::make_unique<Backend>(scale))
{
}

template<class Tout, class Tin, class Backend>
DOp<Tout, Tin, Backend>::DOp() :
    mOp(std::make_unique<Backend>())
{
}

template<class Tout, class Tin, class Backend>
void DOp<Tout, Tin, Backend>::applyImpl(Tout& out, Tin& in)
{
    Profiler::RegionFixture<4> fix("ALegendre::Projector::DOp::applyImpl");

    // // differentiate in place
    // mOp->apply(in, in);

    // // FFT
    // mFft->apply(out, in);
}

template<class Tout, class Tin, class Backend>
void DOp<Tout, Tin, Backend>::applyImpl(Tout& out, const Tin& in)
{
    // // tmp storage
    // std::vector<typename Tin::ScalarType> tmp(in.size());

    // // copy
    // for(std::size_t i = 0; i < in.size(); ++i)
    // {
    //     tmp.data()[i] = in.data()[i];
    // }

    // // tmp view
    // mods_t tmpV(tmp.data(), in.size(), in.dims(), in.pointers(), in.indices());

    // this->applyImpl(out, tmpV);
}

// explicit instantations
// template class DOp<phys_t, mods_t,
//     QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
//     Cpu::DiffOp<mods_t, mods_t, 0, bwd_t>>;
// template class DOp<phys_t, mods_t,
//     QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
//     Cpu::DiffOp<mods_t, mods_t, 0, bwd_t, dealias_m>>;
// template class DOp<phys_t, mods_t,
//     QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
//     Cpu::DiffOp<mods_t, mods_t, 1, bwd_t>>;
// template class DOp<phys_t, mods_t,
//     QuICC::Fft::Fftw::FftOp<phys_t, mods_t>,
//     Cpu::DiffOp<mods_t, mods_t, 1, bwd_t, dealias_m>>;


} // namespace Projector
} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
