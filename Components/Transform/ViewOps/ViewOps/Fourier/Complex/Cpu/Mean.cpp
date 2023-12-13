
#include <iostream>
#include <complex>

#include "Mean.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Util.hpp"
#include "ViewOps/Fourier/Tags.hpp"
#include "ViewOps/Fourier/Complex/Types.hpp"
#include "Profiler/Interface.hpp"

#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#endif

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
namespace Cpu {

template<class Tout, class Tin, class Direction>
MeanOp<Tout, Tin, Direction>::MeanOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, class Direction>
void MeanOp<Tout, Tin, Direction>::applyImpl(Tout& out, const Tin& in)
{
    static_assert(std::is_same_v<typename Tin::LevelType, DCCSC3D::level>,
        "implementation assumes dense, compressed, sparse");

    Profiler::RegionFixture<4> fix("MeanOp::applyImpl");

    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    const auto M = in.lds();

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

            // linear index (:,j,k)
            std::size_t jk = idn*M;
            std::size_t i = 0;
            if (j == 0)
            {
                out.data()[i+jk] = in.data()[i+jk] * fftScaling;
                ++i;
            }
            for (; i < posM; ++i)
            {
                out.data()[i+jk] = 0.0;
            }
            if (j == 0)
            {
                out.data()[i+jk] = in.data()[i+jk] * fftScaling;
                ++i;
            }
            for (; i < M; ++i)
            {
                out.data()[i+jk] = 0.0;
            }
        }
        ++k;
    }
    for (; k < pointers.size()-1 ; ++k)
    {
        for (std::size_t idn = pointers[k]; idn < pointers[k+1]; ++idn)
        {
            // linear index (:,j,k)
            std::size_t jk = idn*M;
            for (std::size_t i = 0; i < M; ++i)
            {
                out.data()[i+jk] = 0.0;
            }
        }
    }


}

// explicit instantations
template class MeanOp<mods_t, mods_t, fwd_t>;
template class MeanOp<mods_t, mods_t, bwd_t>;

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
