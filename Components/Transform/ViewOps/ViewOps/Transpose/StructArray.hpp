/**
 * @file StructArray.hpp
 * @brief Simple fixed size array that is CUDA compatible
 */
#pragma once

// External includes
//

// Project includes
//
#include "View/ViewMacros.hpp"

namespace QuICC {
namespace Transpose {

/// @brief structArray
/// @tparam Scalar
/// @tparam SIZE
template <class Scalar, int SIZE> struct structArray
{
   /// @brief data
   Scalar _data[SIZE];

   /// @brief typedef for pointed data type
   using ScalarType = Scalar;

   /// @brief read only access element in structArray
   /// @param i index
   /// @return reference to element
   QUICC_CUDA_HOSTDEV const Scalar& operator[](std::size_t i) const
   {
      assert(i < SIZE);
      return _data[i];
   }

   /// @brief access element in structArray
   /// @param i index
   /// @return reference to element
   QUICC_CUDA_HOSTDEV Scalar& operator[](std::size_t i)
   {
      assert(i < SIZE);
      return _data[i];
   }

   /// @brief get size of structArray in number of elements
   /// @return SIZE
   QUICC_CUDA_HOSTDEV std::size_t size() const
   {
      return SIZE;
   }

   /// @brief get raw pointer to memory
   /// @return _data
   QUICC_CUDA_HOSTDEV constexpr Scalar* data() const
   {
      return _data;
   }
};

} // namespace Transpose
} // namespace QuICC
