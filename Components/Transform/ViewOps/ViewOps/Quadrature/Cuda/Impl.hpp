/**
 * @file Impl.hpp
 * @brief cuda implementation of Associated Legendre operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Operator/Binary.hpp"
#include "Memory/MemoryResource.hpp"
#include "Memory/Memory.hpp"
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace Quadrature {
/// @brief Cuda backend namespace
namespace Cuda {

/// @brief Derived classes implement Associated Legendre operators.
/// In practice it boils down to a batched matmul with each batch having different sizes.
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Top operator type
/// @tparam Treatment tag to include scaling due to derivative
template<class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
class ImplOp : public Operator::BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top> {
public:
    /// @brief ctor passing memory resource
    /// @param mem memory resource
    ImplOp(std::shared_ptr<Memory::memory_resource> mem);
    /// @brief Default constructor
    ImplOp() = delete;
    /// @brief dtor
    ~ImplOp() = default;
private:
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    /// @param op operator
    void applyImpl(Tout& out, const Tin& in, const Top& op);
    /// @brief Give access to base class
    friend Operator::BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top>;
    /// @brief memory resource
    /// needs shared ptr for memory pools
    /// note, this must call the dtor last
    /// otherwise we cannot dealloc data
    /// \todo consider removing shared ptr and using singleton
    std::shared_ptr<Memory::memory_resource> _mem;
    /// @brief index typedef
    using IndexType = typename Tin::IndexType;
    /// @brief layer index cache
    Memory::MemBlock<IndexType> _layerIndex;
    /// @brief layer width cache
    Memory::MemBlock<IndexType> _layerWidth;
    /// @brief max layer width cache
    std::uint32_t _N;
    /// @brief A (operator) matrix offset cache
    Memory::MemBlock<IndexType> _offSetA;
    /// @brief B (input) matrix offset cache
    Memory::MemBlock<IndexType> _offSetB;
    /// @brief C (output) matrix offset cache
    Memory::MemBlock<IndexType> _offSetC;
};

} // namespace Cuda
} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
