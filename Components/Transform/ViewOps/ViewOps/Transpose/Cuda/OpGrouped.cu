/**
 * @file OpGrouped.cu
 * @brief Grouped Transpose operations on Views
 */

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "OpGrouped.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {


template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Transpose::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(QuICC::Cuda::isDeviceMemory(in.data()));

   if constexpr (std::is_same_v<Perm, p201_t> &&
                 std::is_same_v<typename Tin::AttributesType,
                    View::S1CLCSC3DJIK> &&
                 std::is_same_v<typename Tout::AttributesType,
                    View::DCCSC3DJIK>)
   {
      // dense transpose
      assert(out.size() == in.size());
      const auto I = in.dims()[0];
      const auto J = in.dims()[1];
      const auto K = in.dims()[2];

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 32;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
      numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
      numBlocks.z = 1;

      details::perm<typename Tout::ScalarType, typename Tin::ScalarType, Perm>
         <<<numBlocks, blockSize, sizeof(std::uint32_t) * (2 * I + K)>>>(out,
            in);
      cudaErrChk(cudaPeekAtLastError());
      cudaErrChk(cudaDeviceSynchronize());
   }
   else if constexpr (std::is_same_v<Perm, p120_t> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::DCCSC3DJIK> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::S1CLCSC3DJIK>)
   {
      // dense transpose
      assert(out.size() == in.size());
      const auto I = out.dims()[0];
      const auto J = out.dims()[1];
      const auto K = out.dims()[2];

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 32;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
      numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
      numBlocks.z = 1;

      details::perm<typename Tout::ScalarType, typename Tin::ScalarType, Perm>
         <<<numBlocks, blockSize, sizeof(std::uint32_t) * (2 * I + K)>>>(out,
            in);
      cudaErrChk(cudaPeekAtLastError());
      cudaErrChk(cudaDeviceSynchronize());
   }
   else if constexpr (std::is_same_v<Perm, p201_t> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::DCCSC3D> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::DCCSC3DJIK>)
   {
      // dense transpose
      assert(out.size() <= in.size()); // input might be padded
      assert(out.size() == out.dims()[0] * out.dims()[1] * out.dims()[2]);

      const auto I = in.dims()[0];
      const auto J = in.dims()[1];
      const auto K = in.dims()[2];

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 32;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
      numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
      numBlocks.z = 1;

      details::perm<typename Tout::ScalarType, typename Tin::ScalarType, Perm>
         <<<numBlocks, blockSize>>>(out, in);
      cudaErrChk(cudaPeekAtLastError());
      cudaErrChk(cudaDeviceSynchronize());
   }
   else if constexpr (std::is_same_v<Perm, p120_t> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::DCCSC3DJIK> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::DCCSC3D>)
   {
      // dense transpose
      assert(out.size() >= in.size()); // output might be padded
      assert(out.size() == out.lds() * out.dims()[1] * out.dims()[2]);
      // perm = [1, 2, 0]
      assert(in.dims()[0] == out.dims()[1]);
      assert(in.dims()[1] == out.dims()[2]);
      assert(in.dims()[2] == out.dims()[0]);
      const auto I = in.dims()[0];
      const auto J = in.dims()[1];
      const auto K = in.dims()[2];

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 32;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
      numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
      numBlocks.z = 1;

      details::perm<typename Tout::ScalarType, typename Tin::ScalarType, Perm>
         <<<numBlocks, blockSize>>>(out, in);
      cudaErrChk(cudaPeekAtLastError());
      cudaErrChk(cudaDeviceSynchronize());
   }
   else
   {
      throw std::logic_error("transpose not implemented");
   }
}


// Explicit instantiations
// FT -> AL
template class Op<View::View<double, View::DCCSC3DJIK>,
   View::View<double, View::DCCSC3D>, p201_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3DJIK>,
   View::View<std::complex<double>, View::DCCSC3D>, p201_t>;
// AL -> FT
template class Op<View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3DJIK>, p120_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3D>,
   View::View<std::complex<double>, View::DCCSC3DJIK>, p120_t>;

// AL -> JW
template class Op<View::View<double, View::DCCSC3DJIK>,
   View::View<double, View::S1CLCSC3DJIK>, p201_t>;
template class Op<View::View<std::complex<double>, View::DCCSC3DJIK>,
   View::View<std::complex<double>, View::S1CLCSC3DJIK>, p201_t>;
// JW -> AL
template class Op<View::View<double, View::S1CLCSC3DJIK>,
   View::View<double, View::DCCSC3DJIK>, p120_t>;
template class Op<View::View<std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<std::complex<double>, View::DCCSC3DJIK>, p120_t>;

} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
