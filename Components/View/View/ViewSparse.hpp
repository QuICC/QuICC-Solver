/**
 * @file ViewUnstructured.hpp
 * @brief A view is a n-dimensional dense or sparse tensor reference to a block of memory.
 * It is not the owner of the memory resource.
 */

#pragma once

// System includes
//
#include <cassert>
#include <stdexcept>
#include <vector>
#include <utility>

// Project includes
//
#include "View/ViewBase.hpp"
#include "View/Attributes.hpp"

namespace QuICC {
namespace Memory {

   /// @brief Specialized template for sparse view with attributes
   /// @tparam Scalar element type
   /// @tparam ...Args attributes
   template <class Scalar, class... Args>
   class View<Scalar, Attributes<Args...>>: public ViewBase<Scalar>
   {
    static_assert(hasLevel<Attributes<Args...>>::value, "Attributes has no level");
    static_assert(hasOrder<Attributes<Args...>>::value, "Attributes has no order");

   public:
      /// @brief typedef for pointed data type
      using ScalarType = Scalar;
      /// @brief typedef for indeces and pointers data type
      using IndexType = typename Attributes<Args...>::index;
      /// @brief typedef for aggregated compression level types, i.e. a std::varint of dense_t, compressed_t, etc
      using LevelType = typename Attributes<Args...>::level;
      /// @brief typedef for aggregated logical index types, i.e. a std::varint of i_t, j_t, etc
      using OrderType = typename Attributes<Args...>::order;

   private:
      /// @brief statically defined view rank
      static constexpr std::size_t _rank = std::variant_size_v<LevelType>;

      /// @brief check that attribute sizes matches
      /// @tparam Scalar
      /// @tparam ...Args
      static_assert(_rank == std::variant_size_v<OrderType>, "attributes size mismatch");

      /// @brief logical dimension c-array
      IndexType _dimensions[_rank];

      /// @brief c-array of pointers
      /// ref.: https://arxiv.org/abs/2202.04305
      ViewBase<IndexType> _pointers[_rank];

      /// @brief c-array of indices
      ViewBase<IndexType> _indices[_rank];

   public:
      /// @brief deleted ctor
      View() = delete;
      /// @brief dtor
      virtual ~View() = default;

      /// @brief Generic constructor
      /// @param data
      /// @param dimensions
      /// @param pointers
      /// @param indices
      View(const span<Scalar>& data,
         const std::array<IndexType, _rank> dimensions,
         const std::array<std::vector<IndexType>, _rank>& pointers,
         const std::array<std::vector<IndexType>, _rank>& indices);

      /// @brief Native types constructor
      /// @param data pointer to actual data
      /// @param size size of raw/compressed data
      /// @param dimensions
      /// @param pointers
      /// @param indices
      View(const Scalar* data,
         const std::size_t size,
         const IndexType* dimensions,
         const ViewBase<IndexType>* pointers,
         const ViewBase<IndexType>* indices);

      /// @brief get view rank
      /// @return _rank
      QUICC_CUDA_HOSTDEV static constexpr std::size_t rank() {return _rank;}

      /// @brief get pointer to logical dimensions c-array
      /// @return pointer to _dimensions
      QUICC_CUDA_HOSTDEV const IndexType* dims() const {return _dimensions;}

      /// @brief get pointer to pointers c-array
      /// @return pointer to _pointers
      QUICC_CUDA_HOSTDEV const ViewBase<IndexType>* pointers() const {return _pointers;};

      /// @brief get pointer to indices c-array
      /// @return pointer to _indices
      QUICC_CUDA_HOSTDEV const ViewBase<IndexType>* indices() const {return _indices;};

      /// @brief 1D accessor
      /// @param index logical index
      /// @return element
      ScalarType& operator()(const IndexType index);

      /// @brief 1D accessor
      /// @param index logical index
      /// @return element with const qualifier
      const ScalarType& operator()(const IndexType index) const;

      /// @brief 2D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @return element
      ScalarType& operator()(const IndexType i, const IndexType j);

      /// @brief 2D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @return element with const qualifier
      const ScalarType& operator()(const IndexType i, const IndexType j) const;

      /// @brief 3D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @param k third logical index
      /// @return element
      ScalarType& operator()(const IndexType i, const IndexType j, const IndexType k);

      /// @brief 3D accessor
      /// @param i first logical index
      /// @param j second logical index
      /// @param k third logical index
      /// @return element with const qualifier
      const ScalarType& operator()(const IndexType i, const IndexType j, const IndexType k) const;

   };

   //
   // Definitions
   //

   // Generic ctor
   template <class Scalar, class... Args>
   View<Scalar, Attributes<Args... >>::View(const span<Scalar>& data,
      const std::array<IndexType, _rank> dimensions,
      const std::array<std::vector<IndexType>, _rank>& pointers,
      const std::array<std::vector<IndexType>, _rank>& indices) :
            ViewBase<Scalar>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
         _pointers[i] = pointers[i];
         _indices[i] = indices[i];
      }
   };

   template <class Scalar, class... Args>
   View<Scalar, Attributes<Args... >>::View(const Scalar* data,
      const std::size_t size,
      const IndexType* dimensions,
      const ViewBase<IndexType>* pointers,
      const ViewBase<IndexType>* indices) :
            ViewBase<Scalar>(const_cast<Scalar*>(data), size)
   {

      for (std::size_t i = 0; i < _rank; ++i)
      {
         /// \todo check space consistency
         // assert(QuICC::Cuda::isDeviceMemory(data) ==
         //    QuICC::Cuda::isDeviceMemory(pointers[i]) ||
         //    pointers[i] == nullptr);

         _dimensions[i] = dimensions[i];
         _pointers[i] = pointers[i];
         _indices[i] = indices[i];
      }
   };

   // 1D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType index) const
   {
      static_assert(_rank == 1, "This method is available only to rank 1 Views");
      assert(index < _dimensions[0]);
      /// \todo this can be improved with a binary search
      for (IndexType i = 0; i < _pointers[0][1]; ++i)
      {
         if (_indices[0][i] == index)
         {
            return this->_data[i];
         }
      }
      throw std::logic_error("Trying to refer to an implicit zero.");
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType index)
   {
      return const_cast<Scalar &>(std::as_const(*this).operator()(index));
   }


   // 2D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType i, const IndexType j) const
   {
      static_assert(_rank == 2, "This method is available only to rank 2 Views");
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);
      if constexpr (std::is_same_v<LevelType, DimLevelType<dense_t,compressed_t>>)
      {
         // CSR
         auto ii = i;
         auto jj = j;
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t>>)
         {
            // column major
            std::swap(ii, jj);
         }
         // row major
         // this can be improved with a binary search
         for (IndexType idx = _pointers[1][ii]; idx < _pointers[1][ii+1]; ++idx)
         {
            if (_indices[1][idx] == jj)
            {
               return this->_data[idx];
            }
         }

         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<compressed_t,dense_t>>)
      {
         /// CSC
         throw std::logic_error("Not implemented yet.");
      }
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType i, const IndexType j)
   {
      return const_cast<Scalar &>(std::as_const(*this).operator()(i, j));
   }

   // 3D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<Args... >>::operator()
      (const IndexType i, const IndexType j, const IndexType k) const
   {
      static_assert(_rank == 3, "This method is available only to rank 3 Views");
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);
      assert(k < _dimensions[2]);

      if constexpr (std::is_same_v<LevelType, DimLevelType<compressed_t,dense_t,dense_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            auto n_col = j + k*_dimensions[1];
            for (IndexType idx = _pointers[0][n_col]; idx < _pointers[0][n_col+1]; ++idx)
            {
               if (_indices[0][idx] == i)
               {
                  return this->_data[idx];
               }
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<dense_t, CSC_t, CSC_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major

            // this can be improved with a binary search
            for (IndexType idn = _pointers[1][k]; idn < _pointers[1][k+1]; ++idn)
            {
               if (_indices[1][idn] == j)
               {
                  return this->_data[i+idn*_dimensions[0]];
               }
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<triDense_t,dense_t,triDense_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            // use the knowledge that the column always start from zero and its size depends on the layer
            auto m = _dimensions[0]-k;
            if (i < m)
            {
               std::size_t kmn = 0;
               for(std::size_t ii = 0; ii < k; ++ii)
               {
                  kmn += _dimensions[0]-ii;
               }
               kmn *= _dimensions[1];
               return this->_data[i + j*m + kmn];
            }

         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<triDense_t,dense_t,triCompressed_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            // use the knowledge that the column always start from zero and its size depends on the layer
            std::size_t kmn = 0;
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               auto m = _dimensions[0]-k;
               auto p = _indices[2][idx];
               if ( p == k && i < m)
               {
                  kmn *= _dimensions[1];
                  return this->_data[i + j*m + kmn];
               }
               kmn += _dimensions[0]-p;
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }

   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<Args... >>::operator()
      (const IndexType i, const IndexType j, const IndexType k)
   {
      return const_cast<Scalar &>(std::as_const(*this).operator()(i, j, k));
   }

} // namespace Memory
} // namespace QuICC
