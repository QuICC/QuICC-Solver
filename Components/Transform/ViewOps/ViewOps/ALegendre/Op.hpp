/**
 * @file Op.hpp
 * @brief projector operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Std/Span.hpp"
#include "Operator/Unary.hpp"
#include "Operator/Binary.hpp"
#include "Memory/MemoryResource.hpp"
#include "Memory/Memory.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {

using namespace QuICC::Operator;
using namespace QuICC::Memory;

using QuICC::Patch::std::span;

/// @brief This class implements a ALegendre differentiation and projection
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam Backend type of operator
template<class Tout, class Tin, class Top, class Backend>
class Op : public UnaryBaseOp<Op<Tout, Tin, Top, Backend>, Tout, Tin> {
public:
    /// @brief internal constructor
    /// @param dimensions
    /// @param layers
    /// \todo template might be MP
    Op(span<const typename Top::IndexType> dimensions, span<const typename Top::IndexType> layers,
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

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
