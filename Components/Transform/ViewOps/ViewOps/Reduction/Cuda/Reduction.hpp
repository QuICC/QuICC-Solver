/**
 * @file Reduction.hpp
 * @brief Reduction operations on Views
 * Allows for a reduction operation.
 */
#pragma once

// External includes
//
#include <cstdint>
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Unary.hpp"

namespace QuICC {
/// @brief namespace for Reduction type operations
namespace Reduction {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;
using namespace QuICC::Memory;

/// @brief Reduction operator
/// @tparam Tout n dimensional view
/// @tparam Tin  n-1 dimensional view
/// @tparam Dir axis to perform reduction on
template <class Tout, class Tin, std::uint32_t Dir>
class Op : public UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>
{
public:
   /// @brief ctor passing memory resource
   /// @param mem memory resource
   Op(std::shared_ptr<memory_resource> mem);
   /// @brief Default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

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
   friend UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   /// \todo consider removing shared ptr and using singleton
   std::shared_ptr<memory_resource> _mem;
   /// @brief index typedef
   using IndexType = typename Tin::IndexType;
   /// @brief layer index cache
   MemBlock<IndexType> _layerIndex;
   /// @brief layer width cache
   MemBlock<IndexType> _layerWidth;
   /// @brief max layer width cache
   std::uint32_t _N;
   /// @brief input offset cache
   MemBlock<IndexType> _offSetIn;
   /// @brief output offset cache
   MemBlock<IndexType> _offSetOut;
};

} // namespace Cuda
} // namespace Reduction
} // namespace QuICC
