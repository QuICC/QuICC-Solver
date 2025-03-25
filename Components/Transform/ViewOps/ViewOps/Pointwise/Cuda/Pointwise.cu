#include <cassert>
#include <complex>
#include <cuda/std/complex>

#include "Cuda/CudaUtil.hpp"
#include "Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"

namespace QuICC {
namespace Pointwise {
namespace Cuda {

namespace details {
// naive implementation
template <class Functor, class Tout, class ...Targs>
__global__ void pointwiseKernel(Functor f, Tout out, Targs... args)
{
   const auto M = out.size();
   const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;
   if (m < M)
   {
      out[m] = f(args[m]...);
   }
}

} // namespace details


template <class Functor, class Tout, class ...Targs>
void Op<Functor, Tout, Targs...>::applyImpl(Tout& out, const Targs&... args)
{
   Profiler::RegionFixture<4> fix("Pointwise::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(((QuICC::Cuda::isDeviceMemory(args.data())) && ...));

   // check Tout and Targs.. match Functor op
   using res_t = std::invoke_result_t<Functor, typename Targs::ScalarType...>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(((out.size() == args.size()) && ... ));

   const auto M = out.size();
   dim3 blockSize;
   dim3 numBlocks;
   blockSize.x = 64;
   numBlocks.x = (M + blockSize.x - 1) / blockSize.x;

   details::pointwiseKernel<Functor, Tout, Targs...>
      <<<numBlocks, blockSize>>>(_f, out, args...);
}


// Explicit instantiations
// tests
template class Op<SquareFunctor<double>, View::ViewBase<double>, View::ViewBase<double>>;
template class Op<Abs2Functor<double>, View::ViewBase<double>, View::ViewBase<std::complex<double>>>;
// JW
template class Op<Abs2Functor<double>, View::View<double, View::DCCSC3DJIK>,
   View::View<std::complex<double>, View::DCCSC3DJIK>>;
// Add Mods
template class Op<AddFunctor<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>>;
template class Op<AddFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>>;
template class Op<AddFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>>;
template class Op<AddFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>>;
template class Op<AddFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>>;
template class Op<AddFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;
// Sub Mods
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>>;
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>>;
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>,
   View::View<cuda::std::complex<double>, View::DCCSC3D>>;
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3DJIK>>;
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>,
   View::View<cuda::std::complex<double>, View::S1CLCSC3D>>;
template class Op<SubFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;
// Dot
template class Op<DotFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;
// Cross
template class Op<CrossCompFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;

} // namespace Cuda
} // namespace Pointwise
} // namespace QuICC
