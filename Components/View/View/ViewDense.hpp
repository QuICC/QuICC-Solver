/**
 * @file View.hpp
 * @brief A view is a n-dimensional dense or sparse tensor reference to a block of memory.
 * It is not the owner of the memory resource.
 */

#pragma once

// System includes
//
#include <cassert>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "View/ViewBase.hpp"
#include "View/Attributes.hpp"

namespace QuICC {
namespace View {

   /// @brief Generic template for a dense (structured) view
   /// @tparam Scalar element type
   /// @tparam ...Args
   template <class Scalar, class... Args>
   class ViewDenseBase;

   /// @brief Specialized template for dense view with attributes
   /// @tparam Scalar element type
   /// @tparam ...Args attributes
   template <class Scalar, class... Args>
   class ViewDenseBase<Scalar, Attributes<Args...>>: public ViewBase<Scalar>
   {
    static_assert(hasLevel<Attributes<Args...>>::value, "Attributes has no level");
    static_assert(hasOrder<Attributes<Args...>>::value, "Attributes has no order");

   public:
      /// @brief typedef for pointed data type
      using ScalarType = Scalar;
      /// @brief typedef for indeces and pointers data type
      using IndexType = typename Attributes<Args...>::index;
      /// @brief typedef for aggregated compression level types, i.e. a std::variant of dense_t, compressed_t, etc
      using LevelType = typename Attributes<Args...>::level;
      /// @brief typedef for aggregated logical index types, i.e. a std::variant of i_t, j_t, etc
      using OrderType = typename Attributes<Args...>::order;

   protected:
      /// @brief statically defined view rank 
      static constexpr std::size_t _rank = std::variant_size_v<LevelType>;
      /// @brief check that attribute sizes matches
      /// @tparam Scalar 
      /// @tparam ...Args 
      static_assert(_rank == std::variant_size_v<OrderType>, "attributes size mismatch");

      /// @brief logical dimension
      IndexType _dimensions[_rank];
      
      /// @brief strides in memory
      IndexType _strides[_rank];

   public:
      /// @brief deleted ctor
      ViewDenseBase() = delete;
      /// @brief dtor
      virtual ~ViewDenseBase() = default;
      
      /// @brief import base ctors
      using ViewBase<Scalar>::ViewBase;

      /// @brief get view rank
      /// @return _rank
      QUICC_CUDA_HOSTDEV static constexpr std::size_t rank() {return _rank;}

      /// @brief get pointer to logical dimensions c-array
      /// @return pointer to _dimensions
      QUICC_CUDA_HOSTDEV const IndexType* dims() const {return _dimensions;}
   };

   using dense1D_t = DimLevelType<dense_t>;

   /// @brief 1D dense view
   /// @tparam Scalar element type
   /// @tparam ...Args attributes
   template <class Scalar, class... Args>
   class View<Scalar, Attributes<dense1D_t, Args...>> : public ViewDenseBase<Scalar, Attributes<dense1D_t, Args...>>
   {
   public:
      /// @brief deleted ctor
      View() = delete;
      
      /// @brief dtor
      ~View() = default;

      /// @brief i/// @brief imported from ViewDenseBase
      using typename ViewDenseBase<Scalar, Attributes<dense1D_t, Args...>>::IndexType;
      
      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense1D_t, Args...>>::_rank;
   
      /// @brief /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense1D_t, Args...>>::_dimensions;

      /// @brief 1D fully dense ctor
      /// @param data span ref, useful to accept braced initializer list
      /// @param dimensions std::array of logical dimensions, useful to accept braced initializer list
      View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions);
      
      /// @brief 1D fully dense ctor
      // @param data span ref, useful to accept braced initializer list
      /// @param dimensions pointer to array of logical dimensions
      View(const span<Scalar>& data, const IndexType* dimensions);

      /// @brief 1D accessor
      /// @param index logical index
      /// @return element
      QUICC_CUDA_HOSTDEV Scalar& operator()(const IndexType index);

      /// @brief 1D accessor
      /// @param index logical index
      /// @return element with const attribute
      QUICC_CUDA_HOSTDEV const Scalar& operator()(const IndexType index) const;
   };

   using dense2D_t = DimLevelType<dense_t, dense_t>;

   /// @brief 2D dense view
   /// @tparam Scalar element type
   /// @tparam ...Args attributes
   template <class Scalar, class... Args>
   class View<Scalar, Attributes<dense2D_t, Args...>> : public ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>
   {
   public:
      /// @brief deleted ctor
      View() = delete;
      /// @brief dtor
      ~View() = default;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using typename ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>::IndexType;
      
      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using typename ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>::OrderType;
      
      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>::_rank;
      
      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>::_dimensions;
      
      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense2D_t, Args...>>::_strides;
      
      /// @brief 2D fully dense ctor
      /// @param data span ref, useful to accept braced initializer list
      /// @param dimensions std::array of logical dimensions, useful to accept braced initializer list
      View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions);

      /// @brief 2D fully dense ctor
      /// @param data span ref, useful to accept braced initializer list
      /// @param dimensions pointer to array of logical dimensions
      View(const span<Scalar>& data, const IndexType* dimensions);

      /// @brief 2D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @return element
      QUICC_CUDA_HOSTDEV Scalar& operator()(const IndexType i, const IndexType j);

      /// @brief 2D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @return element with const attribute
      QUICC_CUDA_HOSTDEV const Scalar& operator()(const IndexType i, const IndexType j) const;
   };

   using dense3D_t = DimLevelType<dense_t, dense_t, dense_t>;

   /// View 3D dense
   template <class Scalar, class... Args>
   class View<Scalar, Attributes<dense3D_t, Args...>> : public ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>
   {
   public:
      /// @brief deleted ctor
      View() = delete;
      /// @brief dtor
      ~View() = default;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using typename ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>::IndexType;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using typename ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>::OrderType;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>::_rank;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>::_dimensions;

      /// @brief imported from ViewDenseBase
      /// @tparam Scalar 
      /// @tparam ...Args 
      using ViewDenseBase<Scalar, Attributes<dense3D_t, Args...>>::_strides;

      /// @brief 3D fully dense ctor
      /// @param data span ref, useful to accept braced initializer list
      /// @param dimensions std::array of logical dimensions, useful to accept braced initializer list
      View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions);

      /// @brief 3D fully dense ctor
      /// @param data span ref, useful to accept braced initializer list
      /// @param dimensions pointer to array of logical dimensions
      View(const span<Scalar>& data, const IndexType* dimensions);

      /// @brief 3D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @param k third logical index
      /// @return element
      QUICC_CUDA_HOSTDEV Scalar& operator()(const IndexType i, const IndexType j, const IndexType k);

      /// @brief 3D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @param k third logical index
      /// @return element with const attribute
      QUICC_CUDA_HOSTDEV const Scalar& operator()(const IndexType i, const IndexType j, const IndexType k) const;
   };


   // 1D dense_t ctor
   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense1D_t, Args... >>::View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions) : 
      ViewDenseBase<Scalar, Attributes<dense1D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
   }

   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense1D_t, Args... >>::View(const span<Scalar>& data, const IndexType* dimensions) : 
      ViewDenseBase<Scalar, Attributes<dense1D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
   }

   // 1D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<dense1D_t, Args... >>::operator()(const IndexType index) const
   {
      // static_assert(_rank == 1, "This method is available only to rank 1 Views");
      assert(index < _dimensions[0]);
      return this->_data[index];
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<dense1D_t, Args... >>::operator()(const IndexType index)
   {
      assert(index < _dimensions[0]);
      return this->_data[index];
   }


   // 2D dense_t ctor
   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense2D_t, Args... >>::View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions) :
      ViewDenseBase<Scalar, Attributes<dense2D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
      if constexpr(std::is_same_v<OrderType, LoopOrderType<i_t, j_t>>)
      {
         _strides[0] = 1;
         _strides[1] = _dimensions[0];
      }
      else if constexpr(std::is_same_v<OrderType, LoopOrderType<j_t, i_t>>)
      {
         _strides[1] = 1;
         _strides[0] = _dimensions[1];
      }
   }

   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense2D_t, Args... >>::View(const span<Scalar>& data, const IndexType* dimensions) :
      ViewDenseBase<Scalar, Attributes<dense2D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
      if constexpr(std::is_same_v<OrderType, LoopOrderType<i_t, j_t>>)
      {
         _strides[0] = 1;
         _strides[1] = _dimensions[0];
      }
      else if constexpr(std::is_same_v<OrderType, LoopOrderType<j_t, i_t>>)
      {
         _strides[1] = 1;
         _strides[0] = _dimensions[1];
      }
   }



   // 2D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<dense2D_t, Args... >>::operator()(const IndexType i, const IndexType j) const
   {
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);

      // dense
      return this->_data[i*_strides[0] + j*_strides[1]];
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<dense2D_t, Args... >>::operator()(const IndexType i, const IndexType j)
   {
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);

      // dense
      return this->_data[i*_strides[0] + j*_strides[1]];
   }

   // 3D dense_t ctor
   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense3D_t, Args... >>::View(const span<Scalar>& data, const std::array<IndexType, _rank> dimensions) :
      ViewDenseBase<Scalar, Attributes<dense3D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
      if constexpr(std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
      {
         _strides[0] = 1;
         _strides[1] = _dimensions[0];
         _strides[2] = _dimensions[0]*_dimensions[1];
      }
      else if constexpr(std::is_same_v<OrderType, LoopOrderType<k_t, j_t, i_t>>)
      {

         _strides[2] = 1;
         _strides[1] = _dimensions[2];
         _strides[0] = _dimensions[2]*_dimensions[1];
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }
   }

   template <class Scalar, class... Args>
   View<Scalar, Attributes<dense3D_t, Args... >>::View(const span<Scalar>& data, const IndexType* dimensions) :
      ViewDenseBase<Scalar, Attributes<dense3D_t, Args... >>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
      }
      if constexpr(std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
      {
         _strides[0] = 1;
         _strides[1] = _dimensions[0];
         _strides[2] = _dimensions[0]*_dimensions[1];
      }
      else if constexpr(std::is_same_v<OrderType, LoopOrderType<k_t, j_t, i_t>>)
      {

         _strides[2] = 1;
         _strides[1] = _dimensions[2];
         _strides[0] = _dimensions[2]*_dimensions[1];
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }
   }

   // 3D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<dense3D_t, Args... >>::operator()
      (const IndexType i, const IndexType j, const IndexType k) const
   {
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);
      assert(k < _dimensions[2]);

      return this->_data[i*_strides[0] + j*_strides[1] + k*_strides[2]];
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<dense3D_t, Args... >>::operator()
      (const IndexType i, const IndexType j, const IndexType k)
   {
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);
      assert(k < _dimensions[2]);

      return this->_data[i*_strides[0] + j*_strides[1] + k*_strides[2]];
   }

} // namespace View
} // namespace QuICC
