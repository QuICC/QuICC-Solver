#include <complex>
#include <cuComplex.h>
#include <iostream>

#include "Diff2D.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Util.hpp"
#include "ViewOps/Fourier/Tags.hpp"
#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"


namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
namespace Cuda {

using namespace QuICC::Memory;

/// @brief Compressed sparse layer 3D tensor (Implicit column major)
using mods_t = View<std::complex<double>, DCCSC3D>;

/// @brief thread coarsening factor
constexpr std::size_t tCF = 8;

namespace details
{
    using namespace QuICC::Transform::Fourier::details;

    /// Cuda kernel
    template<std::size_t Ofi, std::size_t Ofj, std::size_t Osi, std::size_t Osj,
        class Direction, std::uint16_t Treatment>
    __global__ void diff2DKernel(mods_t out, const mods_t in, const double scale)
    {
        assert(out.dims()[0] == in.dims()[0]);
        assert(out.dims()[1] == in.dims()[1]);
        assert(out.dims()[2] == in.dims()[2]);

        constexpr bool isFirstComplex = (Ofi+Ofj) % 2;
        constexpr int sgnOfi = 1 - 2*static_cast<int>((Ofi/2) % 2);
        constexpr int sgnOfj = 1 - 2*static_cast<int>((Ofj/2) % 2);
        constexpr int sgnFirst = sgnOfi * sgnOfj;

        const auto M = in.dims()[0];

        // dealias bounds
        std::size_t nDealias = M;
        if constexpr (Treatment & dealias_m)
        {
            nDealias *= dealias::rule;
        }

        // positive / negative coeff bounds
        const auto negM = M / 2;
        const auto posM = negM + M % 2;
        const auto negDealias = nDealias / 2;
        const auto posDealias = negDealias + nDealias % 2;

        double fftScaling = 1.0;
        if constexpr (std::is_same_v<Direction, fwd_t> &&
            !(Treatment & inverse_m))
        {
            fftScaling = 1.0 / static_cast<double>(M);
        }

        if constexpr (std::is_same_v<Direction, fwd_t> &&
            Treatment & inverse_m)
        {
            fftScaling = static_cast<double>(M);
        }

        cuDoubleComplex cF;
        if constexpr (isFirstComplex)
        {
            cF = {0.0, fftScaling};
        }
        else
        {
            cF = {fftScaling, 0.0};
        }

        constexpr bool isSecondComplex = (Osi+Osj) % 2;
        constexpr int sgnOsi = 1 - 2*static_cast<int>((Osi/2) % 2);
        constexpr int sgnOsj = 1 - 2*static_cast<int>((Osj/2) % 2);
        constexpr int sgnSecond = sgnOsi * sgnOsj;
        cuDoubleComplex cS;
        if constexpr (isSecondComplex)
        {
            cS = {0.0, fftScaling};
        }
        else
        {
            cS = {fftScaling, 0.0};
        }

        // Column major
        // Get columns, rows, layers
        auto pointers = in.pointers()[1];
        auto indices = in.indices()[1];
        auto columns = indices.size();

        const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;

        // Column major
        // for (std::size_t m = 0; m < M; ++m)
        if(m < M)
        {
            cuDoubleComplex tmpFm, tmpSm;
            if(m < posM)
            {
                tmpFm = {static_cast<double>(sgnFirst) *
                    fast_pow<Ofi>(static_cast<double>(m)*scale), 0.0};
                tmpSm = {static_cast<double>(sgnSecond) *
                    fast_pow<Osi>(static_cast<double>(m)*scale), 0.0};
            }
            else
            {
                tmpFm = {static_cast<double>(sgnFirst) *
                    fast_pow<Ofi>(-static_cast<double>(M-m)*scale), 0.0};
                tmpSm = {static_cast<double>(sgnSecond) *
                    fast_pow<Osi>(-static_cast<double>(M-m)*scale), 0.0};
            }

            // dealias
            if constexpr (Treatment & dealias_m)
            {
                if(m >= posDealias && m < M - negDealias)
                {
                    tmpFm = {0.0, 0.0};
                    tmpSm = {0.0, 0.0};
                }
            }

            std::size_t k = 0;
            #pragma unroll
            for (std::size_t nnCol = 0; nnCol < tCF; ++nnCol)
            {
                auto nCol = blockIdx.y * tCF + nnCol;
                if (nCol < columns)
                {
                    // \todo could use bisection
                    while(nCol + 1 > pointers[1+k])
                    {
                        ++k;
                    }
                    auto n = indices[nCol];

                    cuDoubleComplex tmpFk = cuCmul({fast_pow<Ofj>(static_cast<double>(k)*scale), 0.0},  cF);
                    cuDoubleComplex tmpSk = cuCmul({fast_pow<Osj>(static_cast<double>(k)*scale), 0.0},  cS);

                    // linear index (m,n,k)
                    auto index = m + nCol*M;

                    auto tmpF = cuCmul(tmpFm, tmpFk);
                    auto tmpS = cuCmul(tmpSm, tmpSk);
                    auto tmpC = cuCadd(tmpF, tmpS);
                    if constexpr (Treatment & inverse_m)
                    {
                        if (k == 0 && n == 0 && (m == 0 || m == posM))
                        {
                            reinterpret_cast<cuDoubleComplex*>(out.data())[index] =
                                cuDoubleComplex{0.0, 0.0};
                        }
                        else
                        {
                            reinterpret_cast<cuDoubleComplex*>(out.data())[index] =
                                cuCdiv(reinterpret_cast<cuDoubleComplex*>(in.data())[index], tmpC);
                        }
                    }
                    else
                    {
                        reinterpret_cast<cuDoubleComplex*>(out.data())[index] =
                            cuCmul(reinterpret_cast<cuDoubleComplex*>(in.data())[index], tmpC);
                    }


                }
            }
        }
    }
}

template<class Tout, class Tin, std::size_t Ofi, std::size_t Ofj,
    std::size_t Osi, std::size_t Osj, class Direction, std::uint16_t Treatment>
Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>::Diff2DOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, std::size_t Ofi, std::size_t Ofj,
    std::size_t Osi, std::size_t Osj, class Direction, std::uint16_t Treatment>
void Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>::applyImpl(Tout& out, const Tin& in)
{
    static_assert(std::is_same_v<typename Tin::LevelType, DCCSC3D::level>,
        "implementation assumes dense, CSC, CSC");

    Profiler::RegionFixture<4> fix("Diff2DOp::applyImpl");

    assert(QuICC::Cuda::isDeviceMemory(out.data()));

    dim3 blockSize;
    blockSize.x = 64;
    blockSize.y = 1;
    blockSize.z = 1;
    dim3 numBlocks;
    numBlocks.x = (in.dims()[0] + blockSize.x - 1) / blockSize.x;
    auto indices = in.indices()[1];
    auto columns = indices.size();
    numBlocks.y = (columns + tCF - 1) / tCF;
    numBlocks.z = 1;
    details::diff2DKernel<Ofi, Ofj, Osi, Osj, Direction, Treatment>
        <<<numBlocks, blockSize>>>(out, in, mScale);
}

// explicit instantations
template class Diff2DOp<mods_t, mods_t, 1, 0, 0, 0, fwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, fwd_t, inverse_m>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 3, 0, 1, 2, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 1, 0, 3, bwd_t>;
template class Diff2DOp<mods_t, mods_t, 2, 0, 0, 2, bwd_t, inverse_m>;

} // namespace Cuda
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
