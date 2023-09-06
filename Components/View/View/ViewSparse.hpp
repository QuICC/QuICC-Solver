/**
 * @file ViewSparse.hpp
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
#include "View/ViewMaps.hpp"
#include "Std/Cuda/Utility.hpp"

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

      /// @brief leading dimensions size
      /// allows for padding of leading dimension in memory
      /// to be used instead of _dim[leading index]
      /// \todo WARNING not all types are implemented/tested
      IndexType _lds;

      /// @brief c-array of pointers
      /// ref.: https://arxiv.org/abs/2202.04305
      ViewBase<IndexType> _pointers[_rank];

      /// @brief c-array of indices
      ViewBase<IndexType> _indices[_rank];


   public:
      /// @brief default ctor for empty view
      View() = default;
      /// @brief dtor
      virtual ~View() = default;

      /// @brief Generic constructor
      /// @param data
      /// @param dimensions
      /// @param pointers
      /// @param indices
      /// @param lds leading dimension size
      View(const span<Scalar>& data,
         const std::array<IndexType, _rank> dimensions,
         const std::array<std::vector<IndexType>, _rank>& pointers,
         const std::array<std::vector<IndexType>, _rank>& indices,
         const IndexType lds = 0);

      /// @brief Native types constructor
      /// @param data pointer to actual data
      /// @param size size of raw/compressed data
      /// @param dimensions
      /// @param pointers
      /// @param indices
      /// @param lds leading dimension size
      View(const Scalar* data,
         const std::size_t size,
         const IndexType* dimensions,
         const ViewBase<IndexType>* pointers,
         const ViewBase<IndexType>* indices,
         const IndexType lds = 0);

      /// @brief get view rank
      /// @return _rank
      QUICC_CUDA_HOSTDEV static constexpr std::size_t rank() {return _rank;}

      /// @brief get pointer to logical dimensions c-array
      /// @return pointer to _dimensions
      QUICC_CUDA_HOSTDEV const IndexType* dims() const {return _dimensions;}

      /// @brief get in memory leading size
      /// @return _lds value
      QUICC_CUDA_HOSTDEV const IndexType lds() const {return _lds;}

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
      const std::array<std::vector<IndexType>, _rank>& indices,
      const IndexType lds) :
            ViewBase<Scalar>(data.data(), data.size())
   {
      for (std::size_t i = 0; i < _rank; ++i)
      {
         _dimensions[i] = dimensions[i];
         _pointers[i] = pointers[i];
         _indices[i] = indices[i];
      }
      // padded leading dimension makes sense only for dense dimensions
      _lds = 0;
      if constexpr (
         std::is_same_v<i_t,   std::variant_alternative_t<0, OrderType>> ||
         std::is_same_v<iInOrder_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<0, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[0];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else if constexpr (std::is_same_v<j_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<1, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[1];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else if constexpr  (std::is_same_v<k_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<2, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[2];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else
      {
         assert(lds == 0);
      }
      assert(_lds <= this->_size);
   };

   template <class Scalar, class... Args>
   View<Scalar, Attributes<Args... >>::View(const Scalar* data,
      const std::size_t size,
      const IndexType* dimensions,
      const ViewBase<IndexType>* pointers,
      const ViewBase<IndexType>* indices,
      const IndexType lds) :
            ViewBase<Scalar>(const_cast<Scalar*>(data), size)
   {
      /// \todo collapse ctors implementations

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

      // padded leading dimension makes sense only for dense dimensions
      _lds = 0;
      if constexpr (
         std::is_same_v<i_t,   std::variant_alternative_t<0, OrderType>> ||
         std::is_same_v<iInOrder_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<0, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[0];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else if constexpr (std::is_same_v<j_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<1, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[1];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else if constexpr  (std::is_same_v<k_t,   std::variant_alternative_t<0, OrderType>>)
      {
         if constexpr (std::is_same_v<dense_t, std::variant_alternative_t<2, LevelType>>)
         {
            if (lds == 0)
            {
               _lds = _dimensions[2];
            }
            else
            {
               _lds = lds;
            }
         }
      }
      else
      {
         assert(lds == 0);
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
            assert(i < this->_size);
            return this->_data[i];
         }
      }
      throw std::logic_error("Trying to refer to an implicit zero.");
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType index)
   {
      return const_cast<Scalar &>(cuda::std::as_const(*this).operator()(index));
   }


   // 2D accessor
   template <class Scalar, class... Args>
   const Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType i, const IndexType j) const
   {
      static_assert(_rank == 2, "This method is available only to rank 2 Views");
      assert(i < _dimensions[0]);
      assert(j < _dimensions[1]);
      if constexpr (std::is_same_v<LevelType, DimLevelType<sparse_t,compressed_t>>)
      {
         // CSR
         static_assert(std::is_same_v<OrderType, LoopOrderType<i_t, j_t>>, "only default mem layout");
         assert(_pointers[1].data() != nullptr);
         assert(_indices[1].data() != nullptr);
         // this can be improved with a binary search
         for (IndexType idx = _pointers[1][i]; idx < _pointers[1][i+1]; ++idx)
         {
            if (_indices[1][idx] == j)
            {
               assert(idx < this->_size);
               return this->_data[idx];
            }
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<compressed_t, sparse_t>>)
      {
         // CSC
         static_assert(std::is_same_v<OrderType, LoopOrderType<i_t, j_t>>, "only default mem layout");
         assert(_pointers[0].data() != nullptr);
         assert(_indices[0].data() != nullptr);
         // this can be improved with a binary search
         for (IndexType idx = _pointers[0][j]; idx < _pointers[0][j+1]; ++idx)
         {
            if (_indices[0][idx] == i)
            {
               assert(idx < this->_size);
               return this->_data[idx];
            }
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else
      {
         throw std::logic_error("Not implemented yet.");
      }
   }

   template <class Scalar, class... Args>
   Scalar& View<Scalar, Attributes<Args... >>::operator()(const IndexType i, const IndexType j)
   {
      return const_cast<Scalar &>(cuda::std::as_const(*this).operator()(i, j));
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

      if constexpr (std::is_same_v<LevelType, DimLevelType<dense_t,dense_t,compressed_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               if (_indices[2][idx] == k)
               {
                  auto kmn = idx*_dimensions[0]*_dimensions[1];
                  assert(i+j*_dimensions[0]+kmn < this->_size);
                  return this->_data[i+j*_dimensions[0]+kmn];
               }
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<dense_t, compressed_t, sparse_t>>)
      {
         if constexpr (
            std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>> ||
            std::is_same_v<OrderType, LoopOrderType<iInOrder_t, j_t, k_t>>)
         {
            // column major

            // map i if needed
            IndexType ii = i;
            if constexpr (std::is_same_v<iInOrder_t, std::variant_alternative_t<0, OrderType>>)
            {
               ii = mapInOrderIndex(i, _dimensions[0], _lds);
            }

            // this can be improved with a binary search
            for (IndexType idn = _pointers[1][k]; idn < _pointers[1][k+1]; ++idn)
            {
               if (_indices[1][idn] == j)
               {
                  assert(ii+idn*_lds < this->_size);
                  return this->_data[ii+idn*_lds];
               }
            }
         }
         else if constexpr (std::is_same_v<OrderType, LoopOrderType<j_t, i_t, k_t>>)
         {
            // JIK

            // this can be improved with a binary search
            IndexType jj = 0;
            for (IndexType idn = _pointers[1][k]; idn < _pointers[1][k+1]; ++idn)
            {
               if (_indices[1][idn] == j)
               {

                  // compute layer shift
                  IndexType kmn = 0;
                  for (IndexType kk = 0; kk < k; ++kk)
                  {

                     kmn += _dimensions[0]* (_pointers[1][kk+1] - _pointers[1][kk]);
                  }

                  IndexType rowSize = _pointers[1][k+1] - _pointers[1][k];
                  IndexType ijk = i*rowSize + jj + kmn;

                  assert(ijk < this->_size);
                  return this->_data[ijk];
               }
               ++jj;
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<triK_t, dense_t, dense_t>>)
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
               assert(i + j*m + kmn < this->_size);
               return this->_data[i + j*m + kmn];
            }

         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<triK_t, dense_t, compressed_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            // use the knowledge that the column always start from zero and its size depends on the layer
            std::size_t kmn = 0;
            const auto n = _dimensions[1];
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               auto m = _dimensions[0]-k;
               auto p = _indices[2][idx];
               if ( p == k && i < m)
               {
                  kmn *= n;
                  assert(i + j*m + kmn < this->_size);
                  return this->_data[i + j*m + kmn];
               }
               // comulative m
               kmn += _dimensions[0]-p;
            }
         }
         else if constexpr (std::is_same_v<OrderType, LoopOrderType<j_t, i_t, k_t>>)
         {
            // i and j are swapped in memory
            std::size_t kmn = 0;
            const auto n = _dimensions[1];
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               auto m = _dimensions[0]-k;
               auto p = _indices[2][idx];
               if ( p == k && i < m)
               {
                  kmn *= n;
                  assert(j + i*n + kmn < this->_size);
                  return this->_data[j + i*n + kmn];
               }
               // comulative m
               kmn += _dimensions[0]-p;
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<dense_t, triK_t, compressed_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            // use the knowledge that the row always start from zero and its size depends on the layer
            std::size_t kmn = 0;
            const auto m = _dimensions[0];
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               auto n = _dimensions[1]-k;
               auto p = _indices[2][idx];
               if ( p == k && j < n)
               {
                  kmn *= m;
                  assert(i + j*m + kmn < this->_size);
                  return this->_data[i + j*m + kmn];
               }
               // comulative n
               kmn += _dimensions[1]-p;
            }
         }
         else if constexpr (std::is_same_v<OrderType, LoopOrderType<j_t, i_t, k_t>>)
         {
            // i and j are swapped in memory
            std::size_t kmn = 0;
            const auto m = _dimensions[0];
            for (IndexType idx = _pointers[2][0]; idx < _pointers[2][1]; ++idx)
            {
               auto n = _dimensions[1]-k;
               auto p = _indices[2][idx];
               if ( p == k && j < n)
               {
                  kmn *= m;
                  assert(j + i*n + kmn < this->_size);
                  return this->_data[j + i*n + kmn];
               }
               // comulative n
               kmn += _dimensions[1]-p;
            }
         }
         else
         {
            throw std::logic_error("Not implemented yet.");
         }
         throw std::logic_error("Trying to refer to an implicit zero.");
      }
      else if constexpr (std::is_same_v<LevelType, DimLevelType<triK_t, compressed_t, sparse_t>>)
      {
         if constexpr (std::is_same_v<OrderType, LoopOrderType<i_t, j_t, k_t>>)
         {
            // column major
            // use the knowledge that the column always start from zero and its size depends on the layer

            const auto M = _dimensions[0];
            const auto K = _dimensions[2];
            assert(_pointers[1].size() == K+1);

            IndexType kmn = 0;
            for (IndexType idx = _pointers[1][k]; idx < _pointers[1][k+1]; ++idx)
            {
               IndexType currColSize = M - k;
               if (j == _indices[1][idx] && i < currColSize)
               {
                  // compute column shift
                  for (IndexType kk = 0; kk < k; ++kk)
                  {
                     IndexType colSize = M - kk;
                     kmn += colSize* (_pointers[1][kk+1] - _pointers[1][kk]);
                  }
                  assert(i + kmn < this->_size);
                  return this->_data[i + kmn];
               }
               kmn += currColSize;
            }
         }
         else if constexpr (std::is_same_v<OrderType, LoopOrderType<j_t, i_t, k_t>>)
         {
            // JIK
            // use the knowledge that the column always start from zero and its size depends on the layer

            const auto M = _dimensions[0];
            const auto K = _dimensions[2];
            assert(_pointers[1].size() == K+1);

            // compute row shift
            IndexType jj = 0;
            for (IndexType idx = _pointers[1][k]; idx < _pointers[1][k+1]; ++idx)
            {
               IndexType currColSize = M - k;
               IndexType currRowSize = _pointers[1][k+1] - _pointers[1][k];
               if (j == _indices[1][idx] && i < currColSize)
               {

                  // compute layer shift
                  IndexType kmn = 0;
                  for (IndexType kk = 0; kk < k; ++kk)
                  {
                     IndexType colSize = M - kk;
                     kmn += colSize* (_pointers[1][kk+1] - _pointers[1][kk]);
                  }

                  assert(i*currRowSize + jj + kmn < this->_size);
                  return this->_data[i*currRowSize + jj + kmn];
               }
               ++jj;
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
      return const_cast<Scalar &>(cuda::std::as_const(*this).operator()(i, j, k));
   }

} // namespace Memory
} // namespace QuICC
