
#include <complex>
#include <cuComplex.h>
#include <iostream>

#include "Mean.hpp"
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
    template<class Direction>
    __global__ void MeanKernel(mods_t out, const mods_t in, const double scale)
    {
        assert(out.dims()[0] == in.dims()[0]);
        assert(out.dims()[1] == in.dims()[1]);
        assert(out.dims()[2] == in.dims()[2]);

        const auto M = in.dims()[0];

        double fftScaling = 1.0;
        if constexpr (std::is_same_v<Direction, fwd_t> )
        {
            fftScaling = 1.0 / static_cast<double>(M);
        }
        cuDoubleComplex cFftScaling = cuDoubleComplex{fftScaling, 0.0};

        // Column major
        // Get columns, rows, layers
        // const auto N = in.dims()[1];
        const auto negM = M / 2;
        const auto posM = negM + M % 2;
        auto pointers = in.pointers()[1];
        auto indices = in.indices()[1];
        auto columns = indices.size();

        constexpr cuDoubleComplex cZero = cuDoubleComplex{0.0, 0.0};

        const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;

        // Column major
        // for (std::size_t m = 0; m < M; ++m)
        if(m < M)
        {
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

                    // linear index (m,n,k)
                    auto index = m + nCol*M;

                    if (k == 0 && n == 0 && (m == 0 || m == posM))
                    {
                        reinterpret_cast<cuDoubleComplex*>(out.data())[index] = cuCmul(
                            reinterpret_cast<cuDoubleComplex*>(in.data())[index], cFftScaling);
                    }
                    else
                    {
                        reinterpret_cast<cuDoubleComplex*>(out.data())[index] = cZero;
                    }

                }
            }
        }
    }
}

template<class Tout, class Tin, class Direction>
MeanOp<Tout, Tin, Direction>::MeanOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, class Direction>
void MeanOp<Tout, Tin, Direction>::applyImpl(Tout& out, const Tin& in)
{
    static_assert(std::is_same_v<typename Tin::LevelType, DCCSC3D::level>,
        "implementation assumes dense, CSC, CSC");

    Profiler::RegionFixture<4> fix("MeanOp::applyImpl");

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
    details::MeanKernel<Direction><<<numBlocks, blockSize>>>(out, in, mScale);
}

// explicit instantations
template class MeanOp<mods_t, mods_t, fwd_t>;
template class MeanOp<mods_t, mods_t, bwd_t>;

} // namespace Cuda
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
