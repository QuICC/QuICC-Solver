/**
 * @file Op.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "ViewOps/Transpose/Cuda/Packing.hpp"
#include "View/View.hpp"
#include "Cuda/CudaUtil.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {


namespace details
{
template<class TDATA>
__global__ void pack(View::ViewBase<TDATA> buffer, const TDATA* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::Attributes<View::dense2D>> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView)
{

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   if (i < I && j < sendCountsView[i])
   {
      buffer[sendBufferDisplsView[i]+j] = *(in + sendDisplsView[i*J+j]);
   }
}

} // namespace details


template <class TDATA>
void pack(View::ViewBase<TDATA> buffer, const TDATA* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::Attributes<View::dense2D>> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView)
{

   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 16;
   blockSize.y = 64;
   blockSize.z = 1;
   numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.z = 1;

   details::pack<TDATA>
      <<<numBlocks, blockSize>>>(buffer, in , sendCountsView, sendDisplsView, sendBufferDisplsView);

}


namespace details
{
template<class TDATA>
__global__ void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::Attributes<View::dense2D>> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView)
{

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   if (i < I && j < recvCountsView[i])
   {
      *(out + recvDisplsView[i*J+j]) = buffer[recvBufferDisplsView[i]+j];
   }
}

} // namespace details


template <class TDATA>
void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::Attributes<View::dense2D>> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView)
{

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 16;
   blockSize.y = 64;
   blockSize.z = 1;
   numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.z = 1;

   details::pack<TDATA>
      <<<numBlocks, blockSize>>>(out, buffer, recvCountsView, recvDisplsView, recvBufferDisplsView);

}



} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
