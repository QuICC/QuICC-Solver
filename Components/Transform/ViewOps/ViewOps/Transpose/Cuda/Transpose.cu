/**
 * @file Transpose.cu
 * @brief Transpose operations on Views
 */

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "Transpose.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;

namespace details
{

template <class Tout, class Tin, class Perm>
__global__ void perm(Tout out, const Tin in)
{
    static_assert(std::is_same_v<Perm, p201_t> &&
                    std::is_same_v<typename Tin::AttributesType, View::DCCSC3D> &&
                    std::is_same_v<typename Tout::AttributesType, View::DCCSC3DJIK>,
      "Not implemented for other types");

    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    const auto I = in.dims()[0];
    const auto J = in.dims()[1];
    const auto K = in.dims()[2];

    if(i < I && j < J) {
        for (std::size_t k = 0; k < K; ++k) {
            std::size_t ijk = i + j*I + k*I*J;
            // plane is row major
            std::size_t jki = j*K + k + i*J*K;
            assert(ijk < in.size());
            assert(jki < out.size());
            out[jki] = in[ijk];
        }
    }
}

} // namespace details


template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Reduction::Cuda::applyImpl");

    assert(QuICC::Cuda::isDeviceMemory(out.data()));
    assert(QuICC::Cuda::isDeviceMemory(in.data()));

    if (std::is_same_v<Perm, p201_t> &&
        std::is_same_v<typename Tin::AttributesType, View::DCCSC3D> &&
        std::is_same_v<typename Tout::AttributesType, View::DCCSC3DJIK>) {
        // dense transpose
        assert(out.size() == in.size());
        assert(out.size() == out.dims()[0]*out.dims()[1]*out.dims()[2]);

        auto I = in.dims()[0];
        auto J = in.dims()[1];
        auto K = in.dims()[2];

        // setup grid
        dim3 blockSize;
        dim3 numBlocks;

        blockSize.x = 32;
        blockSize.y = 32;
        blockSize.z = 1;
        numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
        numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
        numBlocks.z = 1;

        details::perm<Tout, Tin, p201_t><<<numBlocks, blockSize>>>(out, in);

    }
    else {
        throw std::logic_error("transpose not implemented");
    }

}

// Explicit instantiations
template class Op<View::View<double, View::DCCSC3DJIK>, View::View<double, View::DCCSC3D>, p201_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3DJIK>, View::View<std::complex<double>, View::DCCSC3D>, p201_t>;


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
