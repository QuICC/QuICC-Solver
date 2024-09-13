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
               coo[itCoo++] = {
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr)           // k
               };
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
         for (std::size_t i = 0; i < view.dims()[0]; ++i)
         {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
            {
               // 0 1 2 -> i j k && JIK
               coo[itCoo++] = {
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr)           // k
               };
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
               coo[itCoo++] = {
                  static_cast<int>(ptr),          // k
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx])  // j
                  };
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
         for (std::size_t i = 0; i < view.dims()[0]; ++i)
         {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
            {
               // 2 0 1 -> k i j && JIK
               coo[itCoo++] = {
                  static_cast<int>(ptr),          // k
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx])  // j
               };
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
      auto I = view.dims()[0];
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         // std::size_t heightCol = view.dims()[0] - ptr;
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = ptr; i < I; ++i)
            {
               // 0 1 2 -> i j k
               coo[itCoo++] = {
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr)           // k
                  };
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3DJIK> &&
                      std::is_same_v<Perm, p012_t>)
   {
      std::size_t itCoo = 0;
      auto I = view.dims()[0];
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         // std::size_t heightCol = view.dims()[0] - ptr;
         for (std::size_t i = ptr; i < I; ++i)
         {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
            {
               // 0 1 2 -> i j k && JIK
               coo[itCoo++] = {
                  static_cast<int>(i),            // i
                  static_cast<int>(indices[idx]), // j
                  static_cast<int>(ptr)           // k
               };
            }
         }
      }
      return coo;
   }
   else if constexpr (std::is_same_v<typename Tv::AttributesType, S1CLCSC3D> &&
                      std::is_same_v<Perm, p120_t>)
   {
      std::size_t itCoo = 0;
      auto I = view.dims()[0];
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         // std::size_t heightCol = view.dims()[0] - ptr;
         for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
         {
            for (std::size_t i = ptr; i < I; ++i)
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
      auto I = view.dims()[0];
      for (std::size_t ptr = 0; ptr < pointers.size() - 1; ++ptr)
      {
         // std::size_t heightCol = view.dims()[0] - ptr;
         for (std::size_t i = ptr; i < I; ++i)
         {
            for (std::size_t idx = pointers[ptr]; idx < pointers[ptr + 1]; ++idx)
            {
               // 1 2 0 -> j k i && JIK
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
   else
   {
      throw std::logic_error("getCoo not implemented for this type");
   }
   return {};
}

} // namespace View
} // namespace QuICC
