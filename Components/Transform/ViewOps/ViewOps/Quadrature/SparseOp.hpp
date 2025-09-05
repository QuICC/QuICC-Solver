/**
 * @file SparseOp.hpp
 * @brief transform quadrature based sparse operator
 */
#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Binary.hpp"
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "Std/Span.hpp"
#include "ViewOps/Quadrature/ViewBatchedMatmulUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC {
namespace Transform {
namespace Quadrature {

using QuICC::Patch::std::span;

/// @brief This class implements a quadrature based projection
/// with optional differentiation in the third direction
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam Backend type of operator
template <class Tout, class Tin, class Top, class Backend>
class SparseOp : public Operator::UnaryBaseOp<SparseOp<Tout, Tin, Top, Backend>, Tout, Tin>
{
public:
   /// @brief internal constructor
   /// @param mem memory resource for the operator
   SparseOp(std::shared_ptr<QuICC::Memory::memory_resource> mem);
   /// @brief default constructor
   SparseOp() = delete;
   /// @brief dtor
   ~SparseOp() = default;
   /// @brief Allocate operator storage
   /// @param dimensions
   /// @param layers
   void allocOp(span<const typename Top::IndexType> dimensions,
      span<const typename Top::IndexType> layers);

   /// @brief operator accessor
   /// @return operator view
   Top& getOp();

   /// @brief Compress data storage
   void compress(const typename Top::IndexType nnz);

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief action implementation that might modify the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   // void applyImpl(Tout& out, Tin& in);
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<SparseOp<Tout, Tin, Top, Backend>, Tout, Tin>;

   /// @brief pointer to the backend implementation
   std::unique_ptr<Operator::BinaryOp<Tout, Tin, Top>> mImpl;

   ///\todo make private after fixing ctors
public:
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<Memory::memory_resource> _mem;

   /// \todo make single block?
   /// @brief storage for operators
   Memory::MemBlock<typename Top::ScalarType> _opData;
   Memory::MemBlock<typename Top::IndexType> _opIndices;
   Memory::MemBlock<typename Top::IndexType> _opPointers;
   Memory::MemBlock<typename Top::IndexType> _opCsIndices;
   Memory::MemBlock<typename Top::IndexType> _opCsPointers;

   /// @brief View for the operator
   Top _opView;
};


template <class Tout, class Tin, class Top, class Backend>
SparseOp<Tout, Tin, Top, Backend>::SparseOp(std::shared_ptr<Memory::memory_resource> mem) :
    _mem(mem)
{
   // forward memory resource if needed
   if constexpr (std::is_constructible_v<Backend>)
   {
      mImpl = std::make_unique<Backend>();
   }
   else
   {
      mImpl = std::make_unique<Backend>(mem);
   }
}

template <class Tout, class Tin, class Top, class Backend>
void SparseOp<Tout, Tin, Top, Backend>::allocOp(
   span<const typename Top::IndexType> dimensions,
   span<const typename Top::IndexType> layers)
{
   ///\todo move here in details namespace
   auto meta = getOpMeta<Top>(dimensions, layers);

   using namespace QuICC::Memory;
   using IndexType = typename Top::IndexType;

   // Alloc op storage
   _opData = MemBlock<typename Top::ScalarType>(meta.dataSize, _mem.get());
   _opPointers = MemBlock<IndexType>(meta.pointersSize, _mem.get());
   _opIndices = MemBlock<IndexType>(meta.indicesSize, _mem.get());
   _opCsPointers = MemBlock<IndexType>(meta.csPointersSize, _mem.get());
   _opCsIndices = MemBlock<IndexType>(meta.csIndicesSize, _mem.get());

   // Set op view
   using namespace QuICC::View;
   ViewBase<IndexType> pointers[_opView.rank()];
   ViewBase<IndexType> indices[_opView.rank()];
   pointers[meta.idx] =
      ViewBase<IndexType>(_opPointers.data(), _opPointers.size());
   indices[meta.idx] =
      ViewBase<IndexType>(_opIndices.data(), _opIndices.size());
   pointers[meta.csIdx] =
      ViewBase<IndexType>(_opCsPointers.data(), _opCsPointers.size());
   indices[meta.csIdx] =
      ViewBase<IndexType>(_opCsIndices.data(), _opCsIndices.size());
   _opView =
      Top(_opData.data(), meta.dataSize, dimensions.data(), pointers, indices);

   // Adapter for device data
   tempOnHostMemorySpace converterP(pointers[meta.idx], TransferMode::write);
   tempOnHostMemorySpace converterI(indices[meta.idx], TransferMode::write);

   // Set up pointers / indices for operator
   setIndicesAndPointers<Top>(pointers, indices, dimensions, layers);
}

template <class Tout, class Tin, class Top, class Backend>
void SparseOp<Tout, Tin, Top, Backend>::compress(const typename Top::IndexType nnz)
{
   using ScalarType = typename Top::ScalarType;
   using IndexType = typename Top::IndexType;
   using namespace QuICC::Memory;
   using namespace QuICC::View;

   auto _compressedData = MemBlock<ScalarType>(nnz, _mem.get());
   std::copy(_opData.data(), _opData.data() + nnz, _compressedData.data());
   _opData = std::move(_compressedData);

   auto _compressedCsIndices = MemBlock<IndexType>(nnz, _mem.get());
   std::copy(_opCsIndices.data(), _opCsIndices.data() + nnz, _compressedCsIndices.data());
   _opCsIndices = std::move(_compressedCsIndices);

   IndexType csIdx;
   if(_opView.indices()[0].data() == nullptr)
   {
      csIdx = 1;
   }
   else
   {
      csIdx = 0;
   }
   auto indices = _opView.indices();
   ViewBase<IndexType>& csIndices = const_cast<ViewBase<IndexType>*>(indices)[csIdx];
   csIndices = ViewBase<IndexType>(_opCsIndices.data(), _opCsIndices.size());

   // Set op view
   using namespace QuICC::View;
   auto dimensions = _opView.dims();
   auto pointers = _opView.pointers();
   _opView =
      Top(_opData.data(), _opData.size(), dimensions, pointers, indices);
}

template <class Tout, class Tin, class Top, class Backend>
void SparseOp<Tout, Tin, Top, Backend>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Quadrature::Op::applyImpl");

   // Apply backend
   mImpl->apply(out, in, _opView);
}

template <class Tout, class Tin, class Top, class Backend>
Top& SparseOp<Tout, Tin, Top, Backend>::getOp()
{
   return _opView;
}

} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
