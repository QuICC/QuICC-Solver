/**
 * @file Op.Cu
 * @brief Packing for alltoallv and send/recv
 */

// External includes
//
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "View/View.hpp"
#include "ViewOps/Transpose/Cuda/Packing.hpp"

#define QUICC_MAX_PACK_THREADS 512
namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

/// @brief namespace for packing and unpacking details
namespace details {
template <class TDATA, int SIZE>
__global__ void pack(View::ViewBase<TDATA> buffer,
   structArray<const TDATA*, SIZE> in, const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView, const std::int64_t groupSize)
{
   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   const std::size_t i = blockIdx.y * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.x * blockDim.y + threadIdx.y;
   
   // const std::size_t g = blockIdx.z * blockDim.z + threadIdx.z;
   if (i < I)
   {
      int sendCount = sendCountsView[i] / groupSize;
      if (j < sendCount)
      {
         for (int g = 0; g < groupSize; ++g)
         {

            buffer[g * sendCount + sendBufferDisplsView[i] + j] =
               *(in[g] + sendDisplsView[i * J + j]);
         }
      }
   }
}

} // namespace details

template <class TDATA, int SIZE>
void pack(View::ViewBase<TDATA> buffer, structArray<const TDATA*, SIZE> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView, const std::int64_t groupSize)
{

   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];
   // const auto G = in.size();
   
   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 16;
   blockSize.y = 32;
   blockSize.z = 1;
   numBlocks.x = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.y = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.z = 1;

   
   assert(blockSize.x * blockSize.y <= QUICC_MAX_PACK_THREADS);
   __launch_bounds__(QUICC_MAX_PACK_THREADS);
   details::pack<TDATA><<<numBlocks, blockSize>>>(buffer, in, sendCountsView,
      sendDisplsView, sendBufferDisplsView, groupSize);

   cudaErrChk(cudaGetLastError());
   cudaErrChk(cudaDeviceSynchronize());
}


namespace details {
template <class TDATA, int SIZE>
__global__ void unPack(structArray<TDATA*, SIZE> out,
   const View::ViewBase<TDATA> buffer, const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView, const std::int64_t groupSize)
{

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   const std::size_t i = blockIdx.y * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.x * blockDim.y + threadIdx.y;
   // const std::size_t g = blockIdx.z * blockDim.z + threadIdx.z;

   if (i < I)
   {
      int recvCount = recvCountsView[i] / groupSize;
      if (j < recvCount)
      {
         for (int g = 0; g < groupSize; ++g)
         {
            int recvCount = recvCountsView[i] / groupSize;
            *(out[g] + recvDisplsView[i * J + j]) =
               buffer[g * recvCount + recvBufferDisplsView[i] + j];
         }
      }
   }
}

} // namespace details

template <class TDATA, int SIZE>
void unPack(structArray<TDATA*, SIZE> out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView, const std::int64_t groupSize)
{

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];
   // const auto G = out.size();

   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   blockSize.x = 16;
   blockSize.y = 32;
   blockSize.z = 1;
   numBlocks.x = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.y = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.z = 1;

   assert(blockSize.x * blockSize.y <= QUICC_MAX_PACK_THREADS);
   __launch_bounds__(QUICC_MAX_PACK_THREADS);
   details::unPack<TDATA><<<numBlocks, blockSize>>>(out, buffer, recvCountsView,
      recvDisplsView, recvBufferDisplsView, groupSize);

   cudaErrChk(cudaGetLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

// Explicit instantiations
// >>>>
template void pack(View::ViewBase<int> buffer, structArray<const int*, 1> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void pack(View::ViewBase<double> buffer,
   structArray<const double*, 1> in, const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void pack(View::ViewBase<std::complex<double>> buffer,
   structArray<const std::complex<double>*, 1> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<int*, 1> out, const View::ViewBase<int> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<double*, 1> out,
   const View::ViewBase<double> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<std::complex<double>*, 1> out,
   const View::ViewBase<std::complex<double>> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);

template void pack(View::ViewBase<int> buffer, structArray<const int*, 16> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void pack(View::ViewBase<double> buffer,
   structArray<const double*, 16> in, const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void pack(View::ViewBase<std::complex<double>> buffer,
   structArray<const std::complex<double>*, 16> in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<int*, 16> out,
   const View::ViewBase<int> buffer, const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<double*, 16> out,
   const View::ViewBase<double> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);

template void unPack(structArray<std::complex<double>*, 16> out,
   const View::ViewBase<std::complex<double>> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView,
   const std::int64_t groupSize);
// <<<<


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
