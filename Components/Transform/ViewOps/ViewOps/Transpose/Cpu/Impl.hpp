/**
 * @file Impl.hpp
 * @brief Transpose operations on Views serial implementation
 */
#pragma once

// External includes
//

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cpu backends
namespace Cpu {

namespace details {

/// @brief inclusive scan (prefix sum)
/// @tparam T
/// @param vec
template <class T> void pSum(std::vector<T>& vec)
{
   static_assert(std::is_integral_v<T>, "T must be of integral type");
   assert(vec.size() > 0);
   for (std::size_t i = 1; i < vec.size(); ++i)
   {
      vec[i] = vec[i - 1] + vec[i];
   }
}

/// @brief Implementation of permutation [2, 0, 1]
/// for view type input DCCSC3D and output DCCSC3D
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3D>& out,
   const View::View<Tin, View::DCCSC3D>& in)
{
   // dense transpose
   assert(out.size() <= in.size()); // input might be padded
   assert(out.size() == out.dims()[0] * out.dims()[1] * out.dims()[2]);
   // perm = [2, 0, 1]
   assert(in.dims()[0] == out.dims()[2]);
   assert(in.dims()[1] == out.dims()[0]);
   assert(in.dims()[2] == out.dims()[1]);
   auto Ipad = in.lds();
   auto I = in.dims()[0];
   auto J = in.dims()[1];
   auto K = in.dims()[2];
   for (std::size_t k = 0; k < K; ++k)
   {
      for (std::size_t j = 0; j < J; ++j)
      {
         for (std::size_t i = 0; i < I; ++i)
         {
            std::size_t ijk = i + j * Ipad + k * Ipad * J;
            std::size_t jki = j + k * J + i * J * K;
            assert(ijk < in.size());
            assert(jki < out.size());
            out[jki] = in[ijk];
         }
      }
   }
}

/// @brief Implementation of permutation [1, 2, 0]
/// for view type input DCCSC3D and output DCCSC3D
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::DCCSC3D>& out,
   const View::View<Tin, View::DCCSC3D>& in)
{
   // dense transpose
   assert(out.size() >= in.size()); // output might be padded
   assert(out.size() == out.lds() * out.dims()[1] * out.dims()[2]);
   // perm = [1, 2, 0]
   assert(in.dims()[0] == out.dims()[1]);
   assert(in.dims()[1] == out.dims()[2]);
   assert(in.dims()[2] == out.dims()[0]);
   const auto Ipad = out.lds();
   const auto I = out.dims()[0];
   const auto J = out.dims()[1];
   const auto K = out.dims()[2];
   for (std::size_t k = 0; k < K; ++k)
   {
      for (std::size_t j = 0; j < J; ++j)
      {
         std::size_t i = 0;
         for (; i < I; ++i)
         {
            std::size_t ijk = i + j * Ipad + k * Ipad * J;
            std::size_t jki = j + k * J + i * J * K;
            assert(jki < in.size());
            assert(ijk < out.size());
            out[ijk] = in[jki];
         }
         for (; i < Ipad; ++i)
         {
            std::size_t ijk = i + j * Ipad + k * Ipad * J;
            assert(ijk < out.size());
            out[ijk] = 0.0;
         }
      }
   }
}

/// @brief Implementation of permutation [2, 0, 1]
/// for view type input S1CLCSC3D and output DCCSC3D
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm201(View::View<Tout, View::DCCSC3D>& out,
   const View::View<Tin, View::S1CLCSC3D>& in)
{
   // dense transpose
   assert(out.size() == in.size());
   // perm = [2, 0, 1]
   assert(in.dims()[0] == out.dims()[2]);
   assert(in.dims()[1] == out.dims()[0]);
   assert(in.dims()[2] == out.dims()[1]);
   auto I = in.dims()[0];
   auto J = in.dims()[1];
   auto K = in.dims()[2];
   // access S1CLCSC3D
   // cumulative column height is (with ijk) I*k - sum(i)_0^k
   // iSum shifted by 1
   std::vector<std::uint32_t> iSum(K, 0);
   for (std::size_t i = 2; i < K; ++i)
   {
      iSum[i] = iSum[i - 1] + 1;
   }
   pSum(iSum);
   // cumulative row width (with jki)
   // kSum shifted by 1
   std::vector<std::uint32_t> kSum(I, 0);
   kSum[1] = 1;
   for (std::size_t i = 2; i < I; ++i)
   {
      kSum[i] = std::min(kSum[i - 1] + 1, K);
   }
   assert(kSum[I - 1] <= K);
   pSum(kSum);

   for (std::size_t k = 0; k < K; ++k)
   {
      std::size_t Iloc = I - k;
      for (std::size_t j = 0; j < J; ++j)
      {
         for (std::size_t i = k; i < I; ++i)
         {
            std::size_t ijk = i - k + j * Iloc + (k * I - iSum[k]) * J;
            std::size_t jki = j + k * J + kSum[i] * J;
            assert(ijk < in.size());
            assert(jki < out.size());
            out[jki] = in[ijk];
         }
      }
   }
}

/// @brief Implementation of permutation [1, 2, 0]
/// for view type input DCCSC3D and output S1CLCSC3D
/// @tparam Tout
/// @tparam Tin
/// @param out
/// @param in
template <class Tout, class Tin>
void implPerm120(View::View<Tout, View::S1CLCSC3D>& out,
   const View::View<Tin, View::DCCSC3D>& in)
{
   // dense transpose
   assert(out.size() == in.size());
   // perm = [1, 2, 0]
   assert(in.dims()[0] == out.dims()[1]);
   assert(in.dims()[1] == out.dims()[2]);
   assert(in.dims()[2] == out.dims()[0]);
   auto I = out.dims()[0];
   auto J = out.dims()[1];
   auto K = out.dims()[2];
   // access S1CLCSC3D
   // cumulative column height is (with ijk) I*k - sum(i)_0^k
   // iSum shifted by 1
   std::vector<std::uint32_t> iSum(K, 0);
   for (std::size_t i = 2; i < K; ++i)
   {
      iSum[i] = iSum[i - 1] + 1;
   }
   pSum(iSum);
   // cumulative row width (with jki)
   // kSum shifted by 1
   std::vector<std::uint32_t> kSum(I, 0);
   kSum[1] = 1;
   for (std::size_t i = 2; i < I; ++i)
   {
      kSum[i] = std::min(kSum[i - 1] + 1, K);
   }
   assert(kSum[I - 1] <= K);
   pSum(kSum);

   for (std::size_t k = 0; k < K; ++k)
   {
      std::size_t Iloc = I - k;
      for (std::size_t j = 0; j < J; ++j)
      {
         for (std::size_t i = k; i < I; ++i)
         {
            std::size_t ijk = i - k + j * Iloc + (k * I - iSum[k]) * J;
            std::size_t jki = j + k * J + kSum[i] * J;
            assert(jki < in.size());
            assert(ijk < out.size());
            out[ijk] = in[jki];
         }
      }
   }
}

} // namespace details

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
