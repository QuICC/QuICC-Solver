/**
 * @file Impl.hpp
 * @brief Transpose kernel operations on Views
 */
#pragma once

// External includes
//
#include <cassert>
#include <complex>

// Project includes
//
#include "Cuda/CudaUtil.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"

#define QUICC_MAX_TH_NAIVE 2048

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;

namespace details {

template <class Tout, class Tin, class Perm>
__global__ void perm(View::View<Tout, View::DCCSC3DJIK> out,
   const View::View<Tin, View::DCCSC3D> in)
{
   static_assert(std::is_same_v<Perm, p201_t>,
      "Not implemented for other types");

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   const auto Ipad = in.lds();
   const auto I = in.dims()[0];
   const auto J = in.dims()[1];
   const auto K = in.dims()[2];

   if (i < I && j < J)
   {
      for (std::size_t k = 0; k < K; ++k)
      {
         std::size_t ijk = i + j * Ipad + k * Ipad * J;
         // plane is row major
         std::size_t jki = j * K + k + i * J * K;
         assert(ijk < in.size());
         assert(jki < out.size());
         out[jki] = in[ijk];
      }
   }
}

template <class Tout, class Tin, class Perm>
__global__ void perm(View::View<Tout, View::DCCSC3D> out,
   const View::View<Tin, View::DCCSC3DJIK> in)
{
   static_assert(std::is_same_v<Perm, p120_t>,
      "Not implemented for other types");

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   const auto Kpad = out.lds();
   const auto I = in.dims()[0];
   const auto J = in.dims()[1];
   const auto K = in.dims()[2];

   if (i < I && j < J)
   {
      for (std::size_t k = 0; k < K; ++k)
      {
         // plane is row major
         std::size_t ijk = i * J + j + k * I * J;
         std::size_t kij = k + i * Kpad + j * Kpad * I;
         assert(ijk < in.size());
         assert(kij < out.size());
         out[kij] = in[ijk];
      }
   }
}

/// @brief inclusive scan
/// @tparam T
/// @param K
template <class T> __device__ void pSum(T* vec, std::size_t size)
{
   static_assert(std::is_integral_v<T>, "T must be of integral type");
   assert(size > 0);
   for (std::size_t i = 1; i < size; ++i)
   {
      vec[i] = vec[i - 1] + vec[i];
   }
}

template <class Tout, class Tin, class Perm>
__global__ void __launch_bounds__(QUICC_MAX_TH_NAIVE)
   perm(View::View<Tout, View::DCCSC3DJIK> out,
      const View::View<Tin, View::S1CLCSC3DJIK> in)
{
   static_assert(std::is_same_v<Perm, p201_t>,
      "Not implemented for other types");


   const auto I = in.dims()[0];
   const auto J = in.dims()[1];
   const auto K = in.dims()[2];

   // iSum and kSum in shared mem
   extern __shared__ std::uint32_t shared_mem[];
   std::uint32_t* iSum = &shared_mem[0];
   std::uint32_t* kSum = &shared_mem[K];
   std::uint32_t* kLoc = &shared_mem[K + I];

   // access S1CLCSC3D
   // cumulative column height is (with ijk) I*k - sum(i)_0^k

   /// \todo vvv serial chunk to be optimized vvv
   if (threadIdx.x == 0 && threadIdx.y == 0)
   {
      // iSum shifted by 1
      for (std::size_t i = 0; i < K; ++i)
      {
         iSum[i] = 0;
      }
      for (std::size_t i = 2; i < K; ++i)
      {
         iSum[i] = iSum[i - 1] + 1;
      }
      pSum(iSum, K);
      // cumulative row width (with jki)
      // kSum shifted by 1
      kSum[0] = 0;
      // local row width
      kLoc[0] = 1;
      for (std::size_t i = 1; i < I; ++i)
      {
         kSum[i] = min(kSum[i - 1] + 1, K);
         kLoc[i] = min(kLoc[i - 1] + 1, K);
      }
      pSum(kSum, I);
   }

   __syncthreads();

   /// \todo ^^^ serial chunk to be optimized ^^^

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   for (std::size_t k = 0; k < K; ++k)
   {
      // std::size_t Iloc = I - k; // Local column height
      if (i >= k && i < I && j < J)
      {
         std::size_t ijk = (i - k) * J + j + (k * I - iSum[k]) * J;
         std::size_t jki = j * kLoc[i] + k + kSum[i] * J;
         assert(ijk < in.size());
         assert(jki < out.size());
         out[jki] = in[ijk];
      }
   }
}

template <class Tout, class Tin, class Perm>
__global__ void __launch_bounds__(QUICC_MAX_TH_NAIVE)
   perm(View::View<Tout, View::S1CLCSC3DJIK> out,
      const View::View<Tin, View::DCCSC3DJIK> in)
{
   static_assert(std::is_same_v<Perm, p120_t>,
      "Not implemented for other types");


   const auto I = out.dims()[0];
   const auto J = out.dims()[1];
   const auto K = out.dims()[2];

   // iSum and kSum in shared mem
   extern __shared__ std::uint32_t shared_mem[];
   std::uint32_t* iSum = &shared_mem[0];
   std::uint32_t* kSum = &shared_mem[K];
   std::uint32_t* kLoc = &shared_mem[K + I];

   // access S1CLCSC3D
   // cumulative column height is (with ijk) I*k - sum(i)_0^k

   /// \todo vvv serial chunk to be optimized vvv
   if (threadIdx.x == 0 && threadIdx.y == 0)
   {
      // iSum shifted by 1
      for (std::size_t i = 0; i < K; ++i)
      {
         iSum[i] = 0;
      }
      for (std::size_t i = 2; i < K; ++i)
      {
         iSum[i] = iSum[i - 1] + 1;
      }
      pSum(iSum, K);
      // cumulative row width (with jki)
      // kSum shifted by 1
      kSum[0] = 0;
      // local row width
      kLoc[0] = 1;
      for (std::size_t i = 1; i < I; ++i)
      {
         kSum[i] = min(kSum[i - 1] + 1, K);
         kLoc[i] = min(kLoc[i - 1] + 1, K);
      }
      pSum(kSum, I);
   }

   __syncthreads();

   /// \todo ^^^ serial chunk to be optimized ^^^

   const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
   const std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;

   for (std::size_t k = 0; k < K; ++k)
   {
      // std::size_t Iloc = I - k;  // Local column height
      if (i >= k && i < I && j < J)
      {
         std::size_t ijk = (i - k) * J + j + (k * I - iSum[k]) * J;
         std::size_t jki = j * kLoc[i] + k + kSum[i] * J;
         assert(jki < in.size());
         assert(ijk < out.size());
         out[ijk] = in[jki];
      }
   }
}

/// @brief Implementation of permutation [2, 0, 1]
/// for view type input DCCSC3D and output DCCSC3DJIK
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3DJIK>& out,
   const View::View<Tin, View::DCCSC3D>& in);


/// @brief Implementation of permutation [1, 2, 0]
/// for view type input DCCSC3DJIK and output DCCSC3D
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::DCCSC3D>& out,
   const View::View<Tin, View::DCCSC3DJIK>& in);


/// @brief Implementation of permutation [2, 0, 1]
/// for view type input S1CLCSC3DJIK and output DCCSC3DJIK
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3DJIK>& out,
   const View::View<Tin, View::S1CLCSC3DJIK>& in);

/// @brief Implementation of permutation [1, 2, 0]
/// for view type input DCCSC3DJIK and output S1CLCSC3DJIK
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::S1CLCSC3DJIK>& out,
   const View::View<Tin, View::DCCSC3DJIK>& in);


template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3DJIK>& out,
   const View::View<Tin, View::DCCSC3D>& in)
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

   details::perm<Tout, Tin, p201_t><<<numBlocks, blockSize>>>(out, in);
   cudaErrChk(cudaPeekAtLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::DCCSC3D>& out,
   const View::View<Tin, View::DCCSC3DJIK>& in)
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

   details::perm<Tout, Tin, p120_t><<<numBlocks, blockSize>>>(out, in);
   cudaErrChk(cudaPeekAtLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3DJIK>& out,
   const View::View<Tin, View::S1CLCSC3DJIK>& in)

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

   details::perm<Tout, Tin, p201_t>
      <<<numBlocks, blockSize, sizeof(std::uint32_t) * (2 * I + K)>>>(out, in);
   cudaErrChk(cudaPeekAtLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::S1CLCSC3DJIK>& out,
   const View::View<Tin, View::DCCSC3DJIK>& in)
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

   details::perm<Tout, Tin, p120_t>
      <<<numBlocks, blockSize, sizeof(std::uint32_t) * (2 * I + K)>>>(out, in);
   cudaErrChk(cudaPeekAtLastError());
   cudaErrChk(cudaDeviceSynchronize());
}

} // namespace details


} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
