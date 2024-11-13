/**
 * @file Op.hpp
 * @brief Slicewise operations on Views
 * Allows for any user defined Slicewise operation.
 * The operation is defined via a functor object.
 * Value semantic lets a (good) compiler easily inline and
 * remove the indirect call.
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Operator/Nary.hpp"
#include "Memory/MemoryResource.hpp"
#include "Memory/Memory.hpp"

namespace QuICC {
/// @brief namespace for Slicewise type operations
namespace Slicewise {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;

/// @brief Slicewise operator
/// @tparam Functor Nary scalar functor
/// @tparam Tout output View
/// @tparam ...Targs input Views
template <std::uint8_t Dir, class GridBuilder, class Functor, class Tout, class... Targs>
class Op : public NaryBaseOp<Op<Dir, GridBuilder, Functor, Tout, Targs...>, Tout, Targs...>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Op(Functor f, std::shared_ptr<Memory::memory_resource> mem) : _f(f), _mem(mem){};
   /// @brief default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief give access to base class
   friend NaryBaseOp<Op<Dir, GridBuilder, Functor, Tout, Targs...>, Tout, Targs...>;
   /// @brief action implementation
   /// @param out output View
   /// @param ...args input Views
   void applyImpl(Tout& out, const Targs&... args);
   /// @brief specialized implementation for Phi-Theta slice
   void phiThetaImpl(Tout& out, const Targs&... args);
    /// @brief specialized implementation for Phi-R slice
   void phiRImpl(Tout& out, const Targs&... args);
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   /// \todo consider removing shared ptr and using singleton
   std::shared_ptr<Memory::memory_resource> _mem;
   /// @brief index typedef
   using IndexType = typename Tout::IndexType;
   /// @brief layer index cache
   Memory::MemBlock<IndexType> _layerIndex;
   /// @brief layer width cache
   Memory::MemBlock<IndexType> _layerWidth;
   /// @brief max layer width cache
   std::uint32_t _N;
   /// @brief matrix offset cache
   Memory::MemBlock<IndexType> _offSet;
   /// @brief Scalar typedef
   using ScalarType = typename Tout::ScalarType;
   /// @brief reduced grid (only populated layers)
   Memory::MemBlock<ScalarType> _grid;
};

} // namespace Cuda
} // namespace Slicewise
} // namespace QuICC
