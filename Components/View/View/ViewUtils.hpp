/**
 * @file ViewUtils.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <cassert>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace View {

/// @brief compute storage requirements of view
/// @tparam ViewType, either dense or compressed 
/// @param v View
/// @return storage size in bytes
template <class ViewType>
QUICC_CUDA_HOSTDEV std::size_t memBlockSize(ViewType v);

/// @brief compute offset in bytes from data pointer to indices pointer
/// and from data pointer to pointers pointer
/// @tparam ViewType, compressed only
/// @param v View
/// @return array containing the 3 offsets
template <class ViewType>
QUICC_CUDA_HOST std::array<std::size_t, 4> memOffsets(ViewType v);

/// @brief compute offset in bytes from data pointer to indices pointer
/// and from data pointer to pointers pointer assuming they need to be 
/// in the same memblock
/// @tparam ViewType, compressed only
/// @param v View
/// @return array containing the 3 offsets
template <class ViewType>
QUICC_CUDA_HOST std::array<std::size_t, 4> memBlockOffsets(ViewType v);


// 
// Definitions
//

// Storage requirements
template <class ViewType>
std::size_t memBlockSize(ViewType v)
{
   // dense and compressed
   std::size_t size = sizeof(typename ViewType::ScalarType) * v.size();
   
   // compressed only 
   if constexpr (!isLevelTypeFullyDense_v<typename ViewType::LevelType>)
   {
      for (std::size_t i = 0; i < v.rank(); ++i)
      {
            size += v.pointers()[i].size() * sizeof(typename ViewType::IndexType);
            size += v.indices()[i].size() * sizeof(typename ViewType::IndexType);
      }
   }
   return size;
}
   
// offsets
template <class ViewType>
std::array<std::size_t, 4> memOffsets(ViewType v)
{
   static_assert(!isLevelTypeFullyDense_v<typename ViewType::LevelType>);

   std::array<std::size_t, 4> offsets;
   offsets[0] = 0;
   offsets[1] = reinterpret_cast<std::size_t>(v.pointers()) - reinterpret_cast<std::size_t>(v.data());
   offsets[2] = reinterpret_cast<std::size_t>(v.indices()) - reinterpret_cast<std::size_t>(v.data());

   return offsets;
}

// memblock offsets
// no alignement requirements are considered
template <class ViewType>
std::array<std::size_t, 4> memBlockOffsets(ViewType v)
{
   static_assert(!isLevelTypeFullyDense_v<typename ViewType::LevelType>);

   std::array<std::size_t, 4> offsets;
   offsets[0] = 0;
   offsets[1] = sizeof(typename ViewType::ScalarType) * v.size();

   std::size_t size = 0;
   for (std::size_t i = 0; i < v.rank(); ++i)
   {
         size += v.pointers()[i].size() * sizeof(typename ViewType::IndexType);
   }
   offsets[2] = offsets[1] + size;

   return offsets;
}

} // namespace View
} // namespace QuICC
