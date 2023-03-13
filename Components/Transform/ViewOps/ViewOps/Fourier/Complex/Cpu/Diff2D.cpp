
#include <iostream>
#include <complex>

#include "Diff2D.hpp"
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

template<class Tout, class Tin, std::size_t Ofi, std::size_t Ofj,
    std::size_t Osi, std::size_t Osj, class Direction, class Treatment>
Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>::Diff2DOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, std::size_t Ofi, std::size_t Ofj,
    std::size_t Osi, std::size_t Osj, class Direction, class Treatment>
void Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>::applyImpl(Tout& out, const Tin& in)
{
    static_assert(std::is_same_v<typename Tin::LevelType, DCCSC3D::level>,
        "implementation assumes dense, CSC, CSC");

    Profiler::RegionFixture<4> fix("Diff2DOp::applyImpl");

    #ifdef QUICC_USE_CUFFT
    assert(!QuICC::Cuda::isDeviceMemory(out.data()));
    #endif

    using complex_t = typename Tin::ScalarType;
    using float_t = typename Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction,
        Treatment>::ScaleType;

    constexpr bool isFirstComplex = (Ofi+Ofj) % 2;
    constexpr int sgnOfi = 1 - 2*static_cast<int>((Ofi/2) % 2);
    constexpr int sgnOfj = 1 - 2*static_cast<int>((Ofj/2) % 2);
    constexpr int sgnFirst = sgnOfi * sgnOfj;

    const auto M = in.dims()[0];

    float_t fftScaling = 1.0;
    if constexpr (std::is_same_v<Direction, fwd_t> &&
        !std::is_same_v<Treatment, inverse_t>)
    {
        fftScaling = 1.0 / static_cast<float_t>(M);
    }

    if constexpr (std::is_same_v<Direction, fwd_t> &&
        std::is_same_v<Treatment, inverse_t>)
    {
        fftScaling = static_cast<float_t>(M);
    }

    std::conditional_t<isFirstComplex, complex_t, float_t> cF;
    if constexpr (isFirstComplex)
    {
        cF = complex_t(0.0, fftScaling);
    }
    else
    {
        cF = fftScaling;
    }

    constexpr bool isSecondComplex = (Osi+Osj) % 2;
    constexpr int sgnOsi = 1 - 2*static_cast<int>((Osi/2) % 2);
    constexpr int sgnOsj = 1 - 2*static_cast<int>((Osj/2) % 2);
    constexpr int sgnSecond = sgnOsi * sgnOsj;
    std::conditional_t<isSecondComplex, complex_t, float_t> cS;
    if constexpr (isSecondComplex)
    {
        cS = complex_t(0.0, fftScaling);
    }
    else
    {
        cS = fftScaling;
    }

    // Column major
    // Get columns, rows, layers
    const auto negM = M / 2;
    const auto posM = negM + M % 2;
    auto pointers = in.pointers()[1];
    auto indices = in.indices()[1];

    using namespace QuICC::Transform::Fourier::details;

    for (std::size_t k = 0; k < pointers.size()-1 ; ++k)
    {
        for (std::size_t idn = pointers[k]; idn < pointers[k+1]; ++idn)
        {
            std::size_t n = indices[idn];
            // linear index (:,n,k)
            std::size_t nk = idn*M;
            std::size_t m = 0;
            if constexpr (std::is_same_v<Treatment, inverse_t>)
            {
                if (k == 0 && n == 0)
                {
                    out.data()[nk+m] = 0.0;
                    ++m;
                }
            }
            for (; m < posM; ++m)
            {
                auto coeff = (
                    cF * static_cast<float_t>(sgnFirst) *
                        fast_pow<Ofi>(static_cast<float_t>(m)*mScale) *
                        fast_pow<Ofj>(static_cast<float_t>(k)*mScale) +
                    cS * static_cast<float_t>(sgnSecond) *
                        fast_pow<Osi>(static_cast<float_t>(m)*mScale) *
                        fast_pow<Osj>(static_cast<float_t>(k)*mScale));

                if constexpr (std::is_same_v<Treatment, inverse_t>)
                {
                   coeff = 1.0 / coeff;
                }
                out.data()[nk+m] = in.data()[nk+m] * coeff;
            }
            if constexpr (std::is_same_v<Treatment, inverse_t>)
            {
                if (k == 0 && n == 0)
                {
                    out.data()[nk+m] = 0.0;
                    ++m;
                }
            }
            for (; m < M; ++m)
            {
                auto coeff = (
                    cF * static_cast<float_t>(sgnFirst) *
                        fast_pow<Ofi>(-static_cast<float_t>(M-m)*mScale) *
                        fast_pow<Ofj>(static_cast<float_t>(k)*mScale) +
                    cS * static_cast<float_t>(sgnSecond) *
                        fast_pow<Osi>(-static_cast<float_t>(M-m)*mScale) *
                        fast_pow<Osj>(static_cast<float_t>(k)*mScale));

                if constexpr (std::is_same_v<Treatment, inverse_t>)
                {
                   coeff = 1.0 / coeff;
                }
                out.data()[nk+m] = in.data()[nk+m] * coeff;
            }
        }
    }
}

// explicit instantations
using mods_t = View<std::complex<double>, DCCSC3D>;
template class Diff2DOp<mods_t, mods_t, 1, 0, 0, 0, fwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 3, 0, 1, 2, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 1, 0, 3, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_t>;

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
