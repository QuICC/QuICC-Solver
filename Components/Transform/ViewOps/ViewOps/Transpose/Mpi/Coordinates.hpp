/**
 * @file Coordinates.hpp
 * @brief Methods to get absolute coordinate indices
 */
#pragma once

// External includes
//
#include <array>
#include <cassert>
#include <mpi.h>
#include <vector>

// Project includes
//
#include "View/View.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"


namespace QuICC {
namespace View {

/// @brief triplet of indices
using point_t = std::array<int, 3>;

/// @brief Get coordinates (i.e. triplet of indices)
/// for each degree of freedom
/// @tparam Tv
/// @tparam Perm tag describing logical ordering
/// @param view
/// @return
template <class Tv, class Perm> std::vector<point_t> getCoo(const Tv& view)
{
   std::vector<point_t> coo(view.size());
   auto pointers = view.pointers()[1];
   auto indices = view.indices()[1];

   // copy back to cpu
   using namespace QuICC::Memory;
   tempOnHostMemorySpace converterP(pointers, TransferMode::read);
   tempOnHostMemorySpace converterI(indices, TransferMode::read | TransferMode::block);

   using namespace QuICC::Transpose;
   if constexpr (std::is_same_v<typename Tv::AttributesType, DCCSC3D> &&
                 std::is_same_v<Perm, p012_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = 0; i < view.lds(); ++i)
            {
               // 0 1 2 -> i j k
               coo[itCoo++] = {static_cast<int>(i),
                  static_cast<int>(indices[idx]), static_cast<int>(ptr)};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, DCCSC3DJIK> &&
                 std::is_same_v<Perm, p012_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = 0; i < view.lds(); ++i)
            {
               // 0 1 2 -> i j k -> JIK -> j i k
               coo[itCoo++] = {static_cast<int>(indices[idx]), static_cast<int>(i), static_cast<int>(ptr)};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, DCCSC3D> &&
                      std::is_same_v<Perm, p201_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = 0; i < view.lds(); ++i)
            {
               // 2 0 1 -> k i j
               coo[itCoo++] = {static_cast<int>(ptr), static_cast<int>(i),
                  static_cast<int>(indices[idx])};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, DCCSC3DJIK> &&
                      std::is_same_v<Perm, p201_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = 0; i < view.lds(); ++i)
            {
               // 2 0 1 -> k i j -> JIK -> i k j
               coo[itCoo++] = {static_cast<int>(i), static_cast<int>(ptr),
                  static_cast<int>(indices[idx])};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, DCCSC3D> &&
                      std::is_same_v<Perm, p120_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = 0; i < view.lds(); ++i)
            {
               // 1 2 0 -> j k i
               coo[itCoo++] = {
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr),          // k
                  static_cast<int>(i)             // i
               };
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3D> &&
                      std::is_same_v<Perm, p012_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            std::size_t heightCol = view.dims()[0] - ptr;
            for (std::size_t i = 0; i < heightCol; ++i)
            {
               // 0 1 2 -> i j k
               coo[itCoo++] = {static_cast<int>(i),
                  static_cast<int>(indices[idx]), static_cast<int>(ptr)};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3DJIK> &&
                      std::is_same_v<Perm, p012_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            std::size_t heightCol = view.dims()[0] - ptr;
            for (std::size_t i = 0; i < heightCol; ++i)
            {
               // 0 1 2 -> i j k -> JIK -> j i k
               coo[itCoo++] = {static_cast<int>(indices[idx]), static_cast<int>(i), static_cast<int>(ptr)};
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3D> &&
                      std::is_same_v<Perm, p120_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            std::size_t heightCol = view.dims()[0] - ptr;
            for (std::size_t i = 0; i < heightCol; ++i)
            {
               // 1 2 0 -> j k i
               coo[itCoo++] = {
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr),          // k
                  static_cast<int>(i)             // i
               };
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3DJIK> &&
                      std::is_same_v<Perm, p120_t>)
   {
      std::size_t itCoo = 0;
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            std::size_t heightCol = view.dims()[0] - ptr;
            for (std::size_t i = 0; i < heightCol; ++i)
            {
               // 1 2 0 -> j k i -> JIK -> k j i
               coo[itCoo++] = {
                  static_cast<int>(ptr),          // k
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(i)             // i
               };
            }
         }
      }
      return coo;
   }
   else
   {
      throw std::logic_error("getCoo not implemented for this type");
   }
   return {};
}

} // namespace View
} // namespace QuICC
