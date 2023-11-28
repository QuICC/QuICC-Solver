#include <cassert>
#include <complex>

#include "Cuda/CudaUtil.hpp"
#include "Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"

namespace QuICC {
namespace Pointwise {
namespace Cuda {

namespace details {
// naive implementation
template <class Tout, class Tin, class Functor>
__global__ void pointwiseKernel(Tout out, const Tin in, Functor f)
{
   const auto M = out.size();
   const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;
   if (m < M)
   {
      out[m] = f(in[m]);
   }
}

} // namespace details


template <class Tout, class Tin, class Functor>
void Op<Tout, Tin, Functor>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Pointwise::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(QuICC::Cuda::isDeviceMemory(in.data()));

   // check Tout and Tin match Functor op
   using res_t = std::invoke_result_t<Functor, typename Tin::ScalarType>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(out.size() == in.size());

   const auto M = out.size();
   dim3 blockSize;
   dim3 numBlocks;
   blockSize.x = 64;
   numBlocks.x = (M + blockSize.x - 1) / blockSize.x;

   details::pointwiseKernel<Tout, Tin, Functor>
      <<<numBlocks, blockSize>>>(out, in, _f);
}

using namespace QuICC::Memory;

// Explicit instantiations
// tests
template class Op<ViewBase<double>, ViewBase<double>, SquareFunctor<double>>;
template class Op<ViewBase<double>, ViewBase<std::complex<double>>,
   Abs2Functor<double>>;
// JW
template class Op<View<double, DCCSC3DJIK>,
   View<std::complex<double>, DCCSC3DJIK>, Abs2Functor<double>>;

} // namespace Cuda
} // namespace Pointwise
} // namespace QuICC
