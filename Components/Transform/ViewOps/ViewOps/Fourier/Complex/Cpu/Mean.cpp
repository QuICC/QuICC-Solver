
#include <iostream>
#include <complex>

#include "Mean.hpp"
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
namespace Complex {
namespace Cpu {

using namespace QuICC::Memory;

template<class Tout, class Tin, class Direction>
MeanOp<Tout, Tin, Direction>::MeanOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, class Direction>
void MeanOp<Tout, Tin, Direction>::applyImpl(Tout& out, const Tin& in)
{
    static_assert(std::is_same_v<typename Tin::LevelType, DCCSC3D::level>,
        "implementation assumes dense, CSC, CSC");

    Profiler::RegionFixture<4> fix("MeanOp::applyImpl");

    #ifdef QUICC_USE_CUFFT
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    const auto M = in.dims()[0];

    double fftScaling = 1.0;
    if constexpr (std::is_same_v<Direction, fwd_t> )
    {
        fftScaling = 1.0 / static_cast<double>(M);
    }

    // Column major
    // Get columns, rows, layers
    const auto negM = M / 2;
    const auto posM = negM + M % 2;
    auto pointers = in.pointers()[1];
    auto indices = in.indices()[1];

    using namespace QuICC::Transform::Fourier::details;

    std::size_t k = 0;
    {
        for (std::size_t idn = pointers[k]; idn < pointers[k+1]; ++idn)
        {
            std::size_t j = indices[idn];
            std::size_t i = 0;
            if (j == 0)
            {
                out(i,j,k) = in(i,j,k) * fftScaling;
                ++i;
            }
            for (; i < posM; ++i)
            {
                out(i,j,k) = 0.0;
            }
            if (j == 0)
            {
                out(i,j,k) = in(i,j,k) * fftScaling;
                ++i;
            }
            for (; i < M; ++i)
            {
                out(i,j,k) = 0.0;
            }
        }
        ++k;
    }
    for (; k < pointers.size()-1 ; ++k)
    {
        for (std::size_t idn = pointers[k]; idn < pointers[k+1]; ++idn)
        {
            std::size_t j = indices[idn];
            for (std::size_t i = 0; i < M; ++i)
            {
                out(i,j,k) = 0.0;
            }
        }
    }


}

// explicit instantations
using mods_t = View<std::complex<double>, DCCSC3D>;
template class MeanOp<mods_t, mods_t, fwd_t>;
template class MeanOp<mods_t, mods_t, bwd_t>;

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
