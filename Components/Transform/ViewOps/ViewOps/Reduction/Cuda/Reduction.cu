#include <cassert>

#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"
#include "Reduction.hpp"
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC {
namespace Reduction {
namespace Cuda {

using namespace QuICC::View;

namespace details {
// naive implementation of mat reduction
// each thread performs the reduction on a column in the slice
template <class T>
inline __device__ void matReduction(T* out, const T* in, const std::size_t M,
   const std::size_t N)
{
   const std::size_t n = blockIdx.y * blockDim.y + threadIdx.y;

   if (n < N)
   {
      T acc = 0;
      for (std::size_t m = 0; m < M; ++m)
      {
         // in is row major
         auto mn = n + m * N;
         acc += in[mn];
      }
      out[n] = acc;
   }
}

// batched wrapper
// for each z block, extract a slice and call 2D reduction kernel
template <class Tout, class Tin, std::uint32_t Dir>
__global__ void batchedMatReduction(Tout out, const Tin in,
   const ViewBase<std::uint32_t> layerIndex,
   const ViewBase<std::uint32_t> layerWidth,
   const ViewBase<std::uint32_t> offSetOut,
   const ViewBase<std::uint32_t> offSetIn)
{
   static_assert(Dir == 0u &&
                    std::is_same_v<typename Tin::AttributesType, DCCSC3DJIK> &&
                    std::is_same_v<typename Tout::AttributesType, CSC>,
      "Not implemented for other types");

   using IndexType = typename Tin::IndexType;

   const std::size_t l = blockIdx.z;

   // get dimensions
   auto M = in.dims()[0];
   auto N = layerWidth[l];

   // mats pointers
   auto* inPtr = &(in[offSetIn[l]]);
   auto* outPtr = &(out[offSetOut[l]]);

   // compute
   matReduction(outPtr, inPtr, M, N);
}

} // namespace details


template <class Tout, class Tin, std::uint32_t Dir>
Op<Tout, Tin, Dir>::Op(std::shared_ptr<QuICC::Memory::memory_resource> mem) :
    _mem(mem)
{}

template <class Tout, class Tin, std::uint32_t Dir>
void Op<Tout, Tin, Dir>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Reduction::Cuda::applyImpl");

   assert(QuICC::Cuda::isDeviceMemory(out.data()));
   assert(QuICC::Cuda::isDeviceMemory(in.data()));

   using namespace QuICC::Memory;

   // check types consistency
   static_assert(out.rank() == in.rank() - 1, "input/output rank mismatch");
   static_assert(
      std::is_same_v<typename Tin::ScalarType, typename Tout::ScalarType>,
      "input/output scalar type mismatch");

   if constexpr (Dir == 0u &&
                 std::is_same_v<typename Tin::AttributesType, DCCSC3DJIK> &&
                 std::is_same_v<typename Tout::AttributesType, CSC>)
   {
      // check minimal meta data consistency
      assert(out.pointers()[0].size() == in.pointers()[1].size());
      assert(out.indices()[0].size() == in.indices()[1].size());

      using IndexType = typename Tin::IndexType;
      using namespace QuICC::Memory;

      // setup offsets
      if (_layerIndex.data() == nullptr)
      {
         /// \todo move setup to gpu
         auto pointers = in.pointers()[1];
         assert(pointers.data() != nullptr);

         // copy back to cpu for preprocessing
         tempOnHostMemorySpace converterP(pointers,
            TransferMode::read | TransferMode::block);

         _N = 0;
         IndexType nLayers = 0;
         for (IndexType k = 0; k < pointers.size() - 1; ++k)
         {
            IndexType nCols = pointers[k + 1] - pointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
               _N = std::max(_N, nCols);
               ++nLayers;
            }
         }

         // alloc device mem
         assert(_mem.get() != nullptr);
         _layerIndex =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
         _layerWidth =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

         // setup view
         ViewBase<IndexType> vLayerIndex(_layerIndex.data(),
            _layerIndex.size());
         ViewBase<IndexType> vLayerWidth(_layerWidth.data(),
            _layerWidth.size());

         // setup converters
         tempOnHostMemorySpace converterLI(vLayerIndex,
            TransferMode::write | TransferMode::block);
         tempOnHostMemorySpace converterLW(vLayerWidth, TransferMode::write);

         IndexType layCtr = 0;
         for (IndexType k = 0; k < pointers.size() - 1; ++k)
         {
            IndexType nCols = pointers[k + 1] - pointers[k];
            // check if layer is populated
            if (nCols > 0)
            {
               vLayerIndex[layCtr] = k;
               vLayerWidth[layCtr] = nCols;
               ++layCtr;
            }
         }
         // alloc
         _offSetIn =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));
         _offSetOut =
            std::move(QuICC::Memory::MemBlock<IndexType>(nLayers, _mem.get()));

         // setup views
         ViewBase<IndexType> vOffSetIn(_offSetIn.data(), _offSetIn.size());
         ViewBase<IndexType> vOffSetOut(_offSetOut.data(), _offSetOut.size());

         // setup converters
         tempOnHostMemorySpace converterOIn(vOffSetIn, TransferMode::write);
         tempOnHostMemorySpace converterOOut(vOffSetOut, TransferMode::write);

         // exclusive scan offsets
         IndexType M, N;
         vOffSetIn[0] = 0;
         vOffSetOut[0] = 0;
         for (IndexType h = 0; h < nLayers - 1; ++h)
         {
            // get dimensions
            auto M = in.dims()[0];
            auto N = vLayerWidth[h];

            vOffSetIn[h + 1] = vOffSetIn[h] + M * N;
            vOffSetOut[h + 1] = vOffSetOut[h] + N;
         }
      }

      /// \todo more balanced load distribution
      IndexType M = in.dims()[0];
      const IndexType N = _N;
      const IndexType activeLayers = _layerIndex.size();

      // offsets views
      ViewBase<IndexType> layerIndex(_layerIndex.data(), _layerIndex.size());
      ViewBase<IndexType> layerWidth(_layerWidth.data(), _layerWidth.size());
      ViewBase<IndexType> offSetIn(_offSetIn.data(), _offSetIn.size());
      ViewBase<IndexType> offSetOut(_offSetOut.data(), _offSetOut.size());

      // setup grid
      dim3 blockSize;
      dim3 numBlocks;

      blockSize.x = 1;
      blockSize.y = 32;
      blockSize.z = 1;
      numBlocks.x = 1;
      numBlocks.y = (N + blockSize.y - 1) / blockSize.y;
      numBlocks.z = activeLayers;

      details::batchedMatReduction<Tout, Tin, Dir><<<numBlocks, blockSize>>>(
         out, in, layerIndex, layerWidth, offSetOut, offSetIn);
   }
   else
   {
      // reduction not implemented
      static_assert(std::is_same_v<typename Tout::AttributesType, void>,
         "Reduction not implemented for this type");
   }
}

// Explicit instantiations
template class Op<View::View<double, CSC>, View::View<double, DCCSC3DJIK>, 0u>;


} // namespace Cuda
} // namespace Reduction
} // namespace QuICC
