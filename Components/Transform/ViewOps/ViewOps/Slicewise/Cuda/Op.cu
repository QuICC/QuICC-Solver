#include <cassert>
#include <complex>
#include <cuda/std/complex>

#include "Cuda/CudaUtil.hpp"
#include "Slicewise.hpp"
#include "ViewOps/Slicewise/Functors.hpp"

namespace QuICC {
namespace Slicewise {
namespace Cuda {

namespace details {
// naive implementation
template <class Functor, class Tout, class ...Targs>
__global__ void SlicewiseKernel(Functor f, Tout out, Targs... args)
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
   Profiler::RegionFixture<4> fix("Slicewise::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(((QuICC::Cuda::isDeviceMemory(args.data())) && ...));

      // check Tout and Targs.. match Functor op
   using res_t = std::invoke_result_t<Functor, IndexType, typename Targs::ScalarType...>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(((out.size() == args.size()) && ... ));

   // implemented only for physical space
   static_assert(std::is_same_v<Tout, View::View<double, View::DCCSC3D>>);

   // cache populated layers
   auto pointers = out.pointers()[1];
   if (_layerIndex.size() < 1)
   {

      // copy back to cpu for preprocessing
      tempOnHostMemorySpace converterP(pointers, TransferMode::read | TransferMode::block);

      _N = 0;
      IndexType nLayers = 0;
      for (IndexType k = 0; k < pointers.size()-1 ; ++k)
      {
         IndexType nCols = pointers[k+1] - pointers[k];
         // check if layer is populated
         if (nCols > 0)
         {
            _N = std::max(_N, nCols);
            ++nLayers;
         }
      }

      // alloc device mem
      _layerIndex = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
      _layerWidth = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

      // setup view
      ViewBase<IndexType> vLayerIndex(_layerIndex.data(), _layerIndex.size());
      ViewBase<IndexType> vLayerWidth(_layerWidth.data(), _layerWidth.size());

      // setup converters
      tempOnHostMemorySpace converterLI(vLayerIndex, TransferMode::write | TransferMode::block);
      tempOnHostMemorySpace converterLW(vLayerWidth, TransferMode::write);


      IndexType layCtr = 0;
      for (IndexType k = 0; k < modsPointers.size()-1 ; ++k)
      {
         IndexType nCols = modsPointers[k+1] - modsPointers[k];
         // check if layer is populated
         if (nCols > 0)
         {
               vLayerIndex[layCtr] = k;
               vLayerWidth[layCtr] = nCols;
               ++layCtr;
         }
      }

      // alloc
      _offSet = std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

      // setup views
      ViewBase<IndexType> vOffSet(_offSet.data(), _offSet.size());

      // setup converters
      tempOnHostMemorySpace converterO(vOffSet, TransferMode::write);

      // exclusive scan offsets
      vOffSet[0] = 0;
      for (IndexType h = 0; h < nLayers-1 ; ++h)
      {
         vOffSet[h+1] = vOffSet[h] + out.lds()*vLayerWidth[h],;
      }

      // setup grid
      assert(_grid.data() == nullptr);

      // alloc
      _grid = std::move(QuICC::Memory::MemBlock<ScalarType>(nLayers, _mem.get()));

      // setup views
      ViewBase<IndexType> vGrid(_grid.data(), _grid.size());

      // setup converters
      tempOnHostMemorySpace converterO(vGrid, TransferMode::write);


      ::QuICC::Internal::Array igrid;
      ::QuICC::Internal::Array iweights;
      /// \todo template param
      ::QuICC::Polynomial::Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, out.dims()[2]);

      // set grid
      for (IndexType h = 0; h < nLayers-1 ; ++h)
      {
         _grid[h] = igrid[vLayerIndex[h]];
      }
   }

   const auto M = out.size();
   dim3 blockSize;
   dim3 numBlocks;
   blockSize.x = 64;
   numBlocks.x = (M + blockSize.x - 1) / blockSize.x;

   details::SlicewiseKernel<Functor, Tout, Targs...>
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
// Sub Mods
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>,
   View::ViewBase<cuda::std::complex<double>>>;
template class Op<SubFunctor<cuda::std::complex<double>>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>,
   View::View<cuda::std::complex<double>, View::DCCSC3DJIK>>;


} // namespace Cuda
} // namespace Slicewise
} // namespace QuICC
