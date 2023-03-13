#include <complex>
#include <cuComplex.h>
#include <iostream>

#include "Diff.hpp"
#include "View/View.hpp"
#include "ViewOps/Fourier/Util.hpp"
#include "ViewOps/Fourier/Tags.hpp"
#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
namespace Cuda {

using namespace QuICC::Memory;

using mods_t = View<std::complex<double>, DCCSC3D>;

/// @brief thread coarsening factor
constexpr std::size_t tCF = 8;

namespace details
{
    using namespace QuICC::Transform::Fourier::details;

    /// Cuda kernel
    template<std::size_t Order, class Direction, class Treatment>
    __global__ void diffKernel(mods_t out, const mods_t in, const double scale)
    {
        constexpr bool isComplex = Order % 2;
        constexpr int sgn = 1 - 2*static_cast<int>((Order/2) % 2);

        const auto M = out.dims()[0];
        const auto N = out.dims()[1];

        double fftScaling = 1.0;
        if constexpr (std::is_same_v<Direction, fwd_t>)
        {
            fftScaling = 1.0 / static_cast<double>((M-1)*2);
        }

        cuDoubleComplex c;
        if constexpr (isComplex)
        {
            c = {0.0, fftScaling};
        }
        else
        {
            c = {fftScaling, 0.0};
        }

        // map x block and thread to m loop
        const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;

        // Column major
        // Get total number of columns to loop over
        auto indices = in.indices()[1];
        auto columns = indices.size();
        if(m < M)
        {
            cuDoubleComplex tmpR = {static_cast<double>(sgn) *
                fast_pow<Order>(static_cast<double>(m)*scale), 0.0};

            auto tmpC = cuCmul(c, tmpR);

            if constexpr (std::is_same_v<Treatment, zeroP_t>)
            {
                if (m == 0)
                {
                    tmpC = {fftScaling, 0.0};
                }
            }

            if constexpr (std::is_same_v<Treatment, zeroMinusP_t>)
            {
                if (m == 0)
                {
                    tmpC = {-fftScaling, 0.0};
                }
            }

            // map y blocks to columns loop with thread coarsening
            #pragma unroll
            for (std::size_t nn = 0; nn < tCF; ++nn)
            {
                auto n = blockIdx.y * tCF + nn;
                if (n < columns)
                {
                    // linear index (m,n,k)
                    auto index = m + n*M;
                    reinterpret_cast<cuDoubleComplex*>(out.data())[index] =
                        cuCmul(reinterpret_cast<cuDoubleComplex*>(in.data())[index], tmpC);
                }
            }
        }
    }

}

template<class Tout, class Tin, std::size_t Order, class Direction, class Treatment>
DiffOp<Tout, Tin, Order, Direction, Treatment>::DiffOp(ScaleType scale) : mScale(scale){};

template<class Tout, class Tin, std::size_t Order, class Direction, class Treatment>
void DiffOp<Tout, Tin, Order, Direction, Treatment>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("DiffOp::applyImpl");

    assert(out.dims()[0] == in.dims()[0]);
    assert(out.dims()[1] == in.dims()[1]);
    assert(out.dims()[2] == in.dims()[2]);
    assert(QuICC::Cuda::isDeviceMemory(out.data()));

    if constexpr (std::is_same_v<Direction, bwd_t> &&
        std::is_same_v<Treatment, void> &&  Order == 0)
    {
        // if the diff is in place it is a noop
        if(out.data() == in.data())
        {
            return;
        }
    }

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
    details::diffKernel<Order, Direction, Treatment><<<numBlocks, blockSize>>>(out, in, mScale);
}

// explicit instantations
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

} // namespace Cuda
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
