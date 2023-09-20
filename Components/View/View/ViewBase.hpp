/**
 * @file ViewBase.hpp
 * @brief A view is a n-dimensional dense or sparse tensor reference to a block of memory.
 * It is not the owner of the memory resource.
 */

#pragma once

// System includes
//
#include <limits>

// Project includes
//
#include "Std/Span.hpp"
#include "View/ViewMacros.hpp"

using QuICC::Patch::std::span;

namespace QuICC {
/// @brief This namespace provides all the View related code
namespace Memory {

   /// @brief Generic template for a View
   /// @tparam Scalar element type
   /// @tparam ...Args
   template <class Scalar, class... Args>
   class View;

   /** @brief Base class for a view data structure. Represent simply a memory block (pointer and size).
    * Similar to a std::span, it is not the owner of the memory.
    * It is used instead of a std/boost span so that it can be used in device functions.
    * @tparam Scalar element type
    */
   template<class Scalar>
   class ViewBase
   {
   protected:
      /// @brief _data pointer to the memory block
      Scalar* _data{nullptr};
      /// @brief _size number of data elements in the memory block
      std::size_t _size{0};
   public:
      /// @brief typedef for pointed data type
      using ScalarType = Scalar;
      /// @brief ctor, empty ViewBase
      ViewBase() = default;
      /// @brief dtor
      virtual ~ViewBase() = default;
      /// @brief ctor
      /// @param data pointer to memory location
      /// @param size in number of elements
      ViewBase(Scalar* data, std::size_t size) : _data(data), _size(size) {};

      /// @brief copy view from std::vector
      /// @param in input vector
      /// @return *this
      QUICC_CUDA_HOST ViewBase& operator=(const std::vector<Scalar>& in)
      {
         _size = in.size();
         _data = const_cast<Scalar*>(in.data());
         return *this;
      }

      /// @brief access element in ViewBase
      /// @param i index
      /// @return reference to element
      QUICC_CUDA_HOSTDEV  Scalar& operator[](std::size_t i) const
      {
         assert(_data != nullptr);
         assert(i < _size);
         return _data[i];
      }

      /// @brief get size of ViewBase in number of elements
      /// @return _size
      QUICC_CUDA_HOSTDEV std::size_t size() const {return _size;}

      /// @brief get raw pointer to memory
      /// @return _data
      QUICC_CUDA_HOSTDEV constexpr Scalar* data() const {return _data;}
   };

} // namespace Memory
} // namespace QuICC
