/**
 * @file ViewSizeUtils.hpp
 * @brief Utilities to compute data size of given View (type and meta)
 */
#pragma once

// External includes
//
#include <array>
#include <cstdint>
#include <vector>


// Project includes
//
#include "View/Attributes.hpp"
#include "View/ViewUtils.hpp"

namespace QuICC {
namespace View {

/// @brief Get data size of view of type ViewTY
/// from dimensions and meta data
/// @tparam ViewTy
/// @param dim
/// @param meta
/// @return
template <class ViewTy>
std::size_t getDataSize(const std::array<std::uint32_t, 3> dim,
   const ptrAndIdx& meta);

/// @brief Get data size of view of type DCCSC3D
/// from dimensions and meta data
/// @param dim
/// @param meta
/// @return
template <>
std::size_t getDataSize<DCCSC3D>(const std::array<std::uint32_t, 3> dim,
   const ptrAndIdx& meta)
{
   std::size_t cumWidth = 0;
   for (std::size_t ptr = 0; ptr < meta.ptr.size() - 1; ++ptr)
   {
      auto width = meta.ptr[ptr + 1] - meta.ptr[ptr];
      cumWidth += width;
   }
   return cumWidth * dim[0];
}

/// @brief Get data size of view of type S1CLCSC3D
/// from dimensions and meta data
/// @param dim
/// @param meta
/// @return
template <>
std::size_t getDataSize<S1CLCSC3D>(const std::array<std::uint32_t, 3> dim,
   const ptrAndIdx& meta)
{
   std::size_t cumSize = 0;
   for (std::size_t ptr = 0; ptr < meta.ptr.size() - 1; ++ptr)
   {
      auto width = meta.ptr[ptr + 1] - meta.ptr[ptr];
      assert(width >= 0);
      auto height = dim[0] - ptr;
      assert(height > 0);
      cumSize += width * height;
   }
   return cumSize;
}

} // namespace View
} // namespace QuICC
