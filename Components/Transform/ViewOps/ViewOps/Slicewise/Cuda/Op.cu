#include <cassert>
#include <complex>
#include <cuda/std/complex>

#include "Cuda/CudaUtil.hpp"
#include "Op.hpp"
#include "ViewOps/Slicewise/Functors.hpp"
#include "Profiler/Interface.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Casts.hpp"

#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {
namespace Slicewise {
namespace Cuda {

namespace details {
// naive implementation
template <class Functor, class Tout, class ...Targs>
__global__ void SlicewisePhiThetaKernel(
   const View::ViewBase<typename Tout::IndexType> layerWidth,
   const View::ViewBase<typename Tout::IndexType> offSet,
   const View::ViewBase<typename Tout::ScalarType> grid,
   Functor f, Tout out, Targs... args)
{
   const std::size_t l = blockIdx.z;
   const auto M = out.lds();
   const auto N = layerWidth[l];

   const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t n = blockIdx.y * blockDim.y + threadIdx.y;

   if (m < M && n < N)
   {
      auto mnk = offSet[l] + m + n*M;
      out[mnk] = f(grid[l], args[mnk]...);
   }
}

// naive implementation
template <class Functor, class Tout, class ...Targs>
__global__ void SlicewisePhiRKernel(
   const View::ViewBase<typename Tout::ScalarType> grid,
   Functor f, Tout out, Targs... args)
{
   const auto M = out.lds();

   const std::size_t m = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t col = blockIdx.y * blockDim.y + threadIdx.y;

   // column Id index or Theta Idx
   const auto indices = out.indices()[1];
   const auto C = indices.size();

   if (m < M && col < C)
   {
      auto thetaIdx = indices[col];
      auto mnk = m + col*M;
      out[mnk] = f(grid[thetaIdx], args[mnk]...);
   }
}

} // namespace details


template <std::uint8_t Dir, class GridBuilder, class Functor, class Tout, class ...Targs>
void Op<Dir, GridBuilder, Functor, Tout, Targs...>::applyImpl(Tout& out, const Targs&... args)
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
   static_assert(std::is_same_v<Tout, View::View<double, View::DCCSC3D>>,
      "Implemented only for physical space and DCCSC3D layout");

   assert(QuICC::Cuda::isDeviceMemory(out.indices()[1].data()));

   if constexpr(Dir == 1)
   {
      phiRImpl(out, args...);
   }
   else if constexpr(Dir == 2)
   {
      phiThetaImpl(out, args...);
   }
   else
   {
      throw std::logic_error("This slice direction is not implemented.");
   }

}

template <std::uint8_t Dir, class GridBuilder, class Functor, class Tout, class ...Targs>
void Op<Dir, GridBuilder, Functor, Tout, Targs...>::phiRImpl(Tout& out, const Targs&... args)
{
   assert(Dir == 1);

    // setup grid
   if (_grid.data() == nullptr)
   {
      ::QuICC::Internal::Array igrid;
      ::QuICC::Internal::Array iweights;
      GridBuilder quad;
      quad.computeQuadrature(igrid, iweights, out.dims()[Dir]);
      // get theta
      ::QuICC::Internal::Array itheta;
      itheta = igrid.array().acos();

      // alloc
      _grid = std::move(QuICC::Memory::MemBlock<ScalarType>(itheta.size(), _mem.get()));

      // setup views
      View::ViewBase<ScalarType> vGrid(_grid.data(), _grid.size());

      // setup converter
      using namespace QuICC::Memory;
      tempOnHostMemorySpace converterG(vGrid, TransferMode::write);

      // copy
      for (std::size_t i = 0; i < vGrid.size(); ++i)
      {
         vGrid[i] = Internal::cast(itheta[i]);
      }
   }

   // views
   View::ViewBase<ScalarType> grid(_grid.data(), _grid.size());

   const IndexType M = out.size();
   const IndexType N = out.indices()[1].size();

   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 64;
   blockSize.y = 4;
   blockSize.z = 1;
   numBlocks.x = (M + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
   numBlocks.z = 1;

   details::SlicewisePhiRKernel<Functor, Tout, Targs...>
      <<<numBlocks, blockSize>>>(grid, _f, out, args...);

}

template <std::uint8_t Dir, class GridBuilder, class Functor, class Tout, class ...Targs>
void Op<Dir, GridBuilder, Functor, Tout, Targs...>::phiThetaImpl(Tout& out, const Targs&... args)
{
   assert(Dir == 2);

   // cache populated layers
   View::ViewBase<IndexType> pointers = out.pointers()[1];
   if (_layerIndex.size() < 1)
   {

      // copy back to cpu for preprocessing
      using namespace QuICC::Memory;
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
      _layerIndex = std::move(MemBlock<IndexType>(nLayers, _mem.get()));
      _layerWidth = std::move(MemBlock<IndexType>(nLayers, _mem.get()));

      // setup view
      View::ViewBase<IndexType> vLayerIndex(_layerIndex.data(), _layerIndex.size());
      View::ViewBase<IndexType> vLayerWidth(_layerWidth.data(), _layerWidth.size());

      // setup converters
      tempOnHostMemorySpace converterLI(vLayerIndex, TransferMode::write | TransferMode::block);
      tempOnHostMemorySpace converterLW(vLayerWidth, TransferMode::write);


      IndexType layCtr = 0;
      for (IndexType k = 0; k < pointers.size()-1 ; ++k)
      {
         IndexType nCols = pointers[k+1] - pointers[k];
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
      View::ViewBase<IndexType> vOffSet(_offSet.data(), _offSet.size());

      // setup converters
      tempOnHostMemorySpace converterO(vOffSet, TransferMode::write);

      // exclusive scan offsets
      vOffSet[0] = 0;
      for (IndexType h = 0; h < nLayers-1 ; ++h)
      {
         vOffSet[h+1] = vOffSet[h] + out.lds()*vLayerWidth[h];
      }

      // setup grid
      assert(_grid.data() == nullptr);

      // alloc
      _grid = std::move(QuICC::Memory::MemBlock<ScalarType>(nLayers, _mem.get()));

      // setup views
      View::ViewBase<ScalarType> vGrid(_grid.data(), _grid.size());

      // setup converters
      tempOnHostMemorySpace converterG(vGrid, TransferMode::write);

      ::QuICC::Internal::Array igrid;
      ::QuICC::Internal::Array iweights;
      /// \todo template param
      GridBuilder quad;
      quad.computeQuadrature(igrid, iweights, out.dims()[Dir]);

      // set grid
      for (IndexType h = 0; h < nLayers; ++h)
      {
         vGrid[h] = Internal::cast(igrid[vLayerIndex[h]]);
      }
   }

   // views
   View::ViewBase<IndexType> layerWidth(_layerWidth.data(), _layerWidth.size());
   View::ViewBase<IndexType> offSet(_offSet.data(), _offSet.size());
   View::ViewBase<ScalarType> grid(_grid.data(), _grid.size());

   const IndexType M = out.size();
   const IndexType N = _N;
   const IndexType activeLayers = _layerIndex.size();

   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 64;
   blockSize.y = 4;
   blockSize.z = 1;
   numBlocks.x = (M + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
   numBlocks.z = activeLayers;

   details::SlicewisePhiThetaKernel<Functor, Tout, Targs...>
      <<<numBlocks, blockSize>>>(layerWidth, offSet, grid, _f, out, args...);
}

// Explicit instantiations
// physical space
template class Op<2, QuICC::Polynomial::Quadrature::WorlandRule, MulRFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;

template class Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulCosFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;

template class Op<1, QuICC::Polynomial::Quadrature::LegendreRule, MulSinFunctor<double>,
   View::View<double, View::DCCSC3D>,
   View::View<double, View::DCCSC3D>>;

} // namespace Cuda
} // namespace Slicewise
} // namespace QuICC
