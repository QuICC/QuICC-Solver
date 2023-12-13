/**
 * @file Op.hpp
 * @brief transform quadrature based operator
 */
#pragma once

// External includes
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

using namespace QuICC::Operator;
using namespace QuICC::Memory;

using QuICC::Patch::std::span;

/// @brief This class implements a quadrature based projection
/// with optional differentiation in the third direction
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam Backend type of operator
template <class Tout, class Tin, class Top, class Backend>
class Op : public UnaryBaseOp<Op<Tout, Tin, Top, Backend>, Tout, Tin>
{
public:
   /// @brief internal constructor
   /// @param dimensions
   /// @param layers
   /// \todo template might be MP
   Op(span<const typename Top::IndexType> dimensions,
      span<const typename Top::IndexType> layers,
      std::shared_ptr<QuICC::Memory::memory_resource> mem);
   /// @brief default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;
   /// @brief operator accessor
   /// @return operator view
   Top& getOp();

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
   friend UnaryBaseOp<Op<Tout, Tin, Top, Backend>, Tout, Tin>;

   /// @brief pointer to the backend implementation
   std::unique_ptr<BinaryOp<Tout, Tin, Top>> mImpl;

   ///\todo make private after fixing ctors
public:
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<memory_resource> _mem;

   /// \todo make single block?
   /// @brief storage for operators
   MemBlock<typename Top::ScalarType> _opData;
   MemBlock<typename Top::IndexType> _opIndices;
   MemBlock<typename Top::IndexType> _opPointers;

   /// @brief View for the operator
   Top _opView;
};


template <class Tout, class Tin, class Top, class Backend>
Op<Tout, Tin, Top, Backend>::Op(span<const typename Top::IndexType> dimensions,
   span<const typename Top::IndexType> layers,
   std::shared_ptr<QuICC::Memory::memory_resource> mem) :
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

   ///\todo move here in details namespace
   auto meta = getOpMeta<Top>(dimensions, layers);

   using namespace QuICC::Memory;
   using IndexType = typename Top::IndexType;

   // Alloc op storage
   _opData = MemBlock<typename Top::ScalarType>(meta.dataSize, _mem.get());
   _opPointers = MemBlock<IndexType>(meta.pointersSize, _mem.get());
   _opIndices = MemBlock<IndexType>(meta.indicesSize, _mem.get());

   // Set op view
   ViewBase<IndexType> pointers[_opView.rank()];
   ViewBase<IndexType> indices[_opView.rank()];
   pointers[meta.idx] =
      ViewBase<IndexType>(_opPointers.data(), _opPointers.size());
   indices[meta.idx] =
      ViewBase<IndexType>(_opIndices.data(), _opIndices.size());
   _opView =
      Top(_opData.data(), meta.dataSize, dimensions.data(), pointers, indices);

   // Adapter for device data
   tempOnHostMemorySpace converterP(pointers[meta.idx], TransferMode::write);
   tempOnHostMemorySpace converterI(indices[meta.idx], TransferMode::write);

   // Set up pointers / indices for operator
   setIndicesAndPointers<Top>(pointers, indices, dimensions, layers);
}

template <class Tout, class Tin, class Top, class Backend>
void Op<Tout, Tin, Top, Backend>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Quadrature::Projector::Op::applyImpl");

   // Apply backend
   mImpl->apply(out, in, _opView);
}

template <class Tout, class Tin, class Top, class Backend>
Top& Op<Tout, Tin, Top, Backend>::getOp()
{
   return _opView;
}

} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
