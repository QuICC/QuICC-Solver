/**
 * @file Op.Cu
 * @brief Packing for alltoallv and send/recv
 */

// External includes
//
#include <complex>

// Project includes
//
#include "ViewOps/Transpose/Cuda/Packing.hpp"
#include "View/View.hpp"
#include "Cuda/CudaUtil.hpp"

#define SERIAL_PACKING 0

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
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView)
{
   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   #if SERIAL_PACKING
   for (size_t i = 0; i < I; ++i)
   {
      for (int j = 0; j < sendCountsView[i]; ++j)
      {
         buffer[sendBufferDisplsView[i]+j] = *(in + sendDisplsView[i*J+j]);
      }
   }
   #else
   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   if (i < I && j < sendCountsView[i])
   {
      buffer[sendBufferDisplsView[i]+j] = *(in + sendDisplsView[i*J+j]);
   }
   #endif



}

} // namespace details


template <class TDATA>
void pack(View::ViewBase<TDATA> buffer, const TDATA* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView)
{

   const auto I = sendDisplsView.dims()[0];
   const auto J = sendDisplsView.dims()[1];

   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   #if SERIAL_PACKING
   blockSize.x = 1;
   blockSize.y = 1;
   blockSize.z = 1;
   numBlocks.x = 1;
   numBlocks.y = 1;
   numBlocks.z = 1;
   #else
   blockSize.x = 16;
   blockSize.y = 64;
   blockSize.z = 1;
   numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.z = 1;
   #endif

   details::pack<TDATA>
      <<<numBlocks, blockSize>>>(buffer, in, sendCountsView, sendDisplsView, sendBufferDisplsView);

   cudaErrChk(cudaDeviceSynchronize());
}


namespace details
{
template<class TDATA>
__global__ void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView)
{
   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   #if SERIAL_PACKING
   for (size_t i = 0; i < I; ++i)
   {
      for (int j = 0; j < recvCountsView[i]; ++j)
      {
         *(out + recvDisplsView[i*J+j]) = buffer[recvBufferDisplsView[i]+j];
      }
   }
   #else
   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   if (i < I && j < recvCountsView[i])
   {
      *(out + recvDisplsView[i*J+j]) = buffer[recvBufferDisplsView[i]+j];
   }
   #endif

}

} // namespace details


template <class TDATA>
void unPack(TDATA* out, const View::ViewBase<TDATA> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView)
{

   const auto I = recvDisplsView.dims()[0];
   const auto J = recvDisplsView.dims()[1];

   // setup grid
   dim3 blockSize;
   dim3 numBlocks;

   #if SERIAL_PACKING
   blockSize.x = 1;
   blockSize.y = 1;
   blockSize.z = 1;
   numBlocks.x = 1;
   numBlocks.y = 1;
   numBlocks.z = 1;
   #else
   blockSize.x = 16;
   blockSize.y = 64;
   blockSize.z = 1;
   numBlocks.x = (I + blockSize.x - 1) / blockSize.x;
   numBlocks.y = (J + blockSize.y - 1) / blockSize.y;
   numBlocks.z = 1;
   #endif


   details::unPack<TDATA>
      <<<numBlocks, blockSize>>>(out, buffer, recvCountsView, recvDisplsView, recvBufferDisplsView);

   cudaErrChk(cudaDeviceSynchronize());
}

// Explicit instantiations

template
void pack(View::ViewBase<int> buffer, const int* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);

template
void pack(View::ViewBase<double> buffer, const double* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);

template
void pack(View::ViewBase<std::complex<double>> buffer, const std::complex<double>* in,
   const View::ViewBase<int> sendCountsView,
   const View::View<int, View::dense2DRM> sendDisplsView,
   const View::ViewBase<int> sendBufferDisplsView);

template
void unPack(int* out, const View::ViewBase<int> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView);

template
void unPack(double* out, const View::ViewBase<double> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView);

template
void unPack(std::complex<double>* out, const View::ViewBase<std::complex<double>> buffer,
   const View::ViewBase<int> recvCountsView,
   const View::View<int, View::dense2DRM> recvDisplsView,
   const View::ViewBase<int> recvBufferDisplsView);


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
