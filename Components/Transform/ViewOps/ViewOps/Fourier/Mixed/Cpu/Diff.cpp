#include <iostream>
#include <complex>

#include "Diff.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Util.hpp"
#include "ViewOps/Fourier/Tags.hpp"
#include "Profiler/Interface.hpp"

#ifdef QUICC_USE_CUFFT
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
namespace Cpu {

using namespace QuICC::Memory;

template<class Tout, class Tin, std::size_t Order, class Direction, class Treatment>
DiffOp<Tout, Tin, Order, Direction, Treatment>::DiffOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, std::size_t Order, class Direction, class Treatment>
void DiffOp<Tout, Tin, Order, Direction, Treatment>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("DiffOp::applyImpl");

    assert(out.dims()[0] == in.dims()[0]);
    assert(out.dims()[1] == in.dims()[1]);
    assert(out.dims()[2] == in.dims()[2]);

    #ifdef QUICC_USE_CUFFT
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    if constexpr (std::is_same_v<Direction, bwd_t> &&
        std::is_same_v<Treatment, void> &&  Order == 0)
    {
        // if the diff is in place it is a noop
        if(out.data() == in.data())
        {
            return;
        }
    }

    constexpr bool isComplex = Order % 2;
    constexpr int sgn = 1 - 2*static_cast<int>((Order/2) % 2);

    using complex_t = typename Tin::ScalarType;
    using float_t = typename DiffOp<Tout, Tin, Order, Direction, Treatment>::ScaleType;
    std::conditional_t<isComplex, complex_t, float_t> c;

    const auto M = in.dims()[0];

    float_t fftScaling = 1.0;
    if constexpr (std::is_same_v<Direction, fwd_t>)
    {
        fftScaling = 1.0 / static_cast<float_t>((M-1)*2);
    }

    if constexpr (isComplex)
    {
        c = complex_t(0.0, fftScaling);
    }
    else
    {
        c = fftScaling;
    }

    using namespace QuICC::Transform::Fourier::details;

    // Column major
    // Get total number of columns to loop over
    auto indices = in.indices()[1];
    auto columns = indices.size();
    for (std::size_t col = 0; col < columns ; ++col)
    {
        // linear index (:,n,k)
        std::size_t nk = M*col;
        std::size_t m = 0;
        if constexpr (std::is_same_v<Treatment, zeroP_t>)
        {
            out.data()[nk+m] = in.data()[nk+m] * fftScaling;
            ++m;
        }

        if constexpr (std::is_same_v<Treatment, zeroMinusP_t>)
        {
            out.data()[nk+m] = -in.data()[nk+m] * fftScaling;
            ++m;
        }

        for (; m < M; ++m)
        {
            out.data()[nk+m] = in.data()[nk+m] * c * static_cast<float_t>(sgn) *
                fast_pow<Order>(static_cast<float_t>(m)*mScale);
        }
    }

}

// explicit instantations
using mods_t = View<std::complex<double>, DCCSC3D>;
template class DiffOp<mods_t, mods_t, 0, fwd_t>;
template class DiffOp<mods_t, mods_t, 1, fwd_t>;
template class DiffOp<mods_t, mods_t, 1, fwd_t, zeroP_t>;
template class DiffOp<mods_t, mods_t, 1, fwd_t, zeroMinusP_t>;
template class DiffOp<mods_t, mods_t, 2, fwd_t>;
template class DiffOp<mods_t, mods_t, 3, fwd_t>;
template class DiffOp<mods_t, mods_t, 4, fwd_t>;
template class DiffOp<mods_t, mods_t, 0, bwd_t>;
template class DiffOp<mods_t, mods_t, 1, bwd_t>;
template class DiffOp<mods_t, mods_t, 2, bwd_t>;
template class DiffOp<mods_t, mods_t, 3, bwd_t>;
template class DiffOp<mods_t, mods_t, 4, bwd_t>;

} // namespace Cpu
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
