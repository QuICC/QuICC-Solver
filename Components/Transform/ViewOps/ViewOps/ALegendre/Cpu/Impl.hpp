/**
 * @file Impl.hpp
 * @brief cpu implementation of Associated Legendre operator
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
namespace ALegendre {
/// @brief Cpu backend namespace
namespace Cpu {


using namespace QuICC::Operator;
using namespace QuICC::Memory;

/// @brief Derived classes implement Associate Legendre operators.
/// In practice it boils down to a batched matmul with each batch having different sizes.
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Top operator type
/// @tparam Treatment tag to include scaling due to derivative
template<class Tout, class Tin, class Top, std::uint16_t Treatment = 0>
class ImplOp : public BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top> {
public:
    /// @brief Default constructor
    ImplOp() = default;
    /// @brief dtor
    ~ImplOp() = default;
private:
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    /// @param op operator
    void applyImpl(Tout& out, const Tin& in, const Top& op);
    /// @brief Give access to base class
    friend BinaryBaseOp<ImplOp<Tout, Tin, Top, Treatment>, Tout, Tin, Top>;
    /// @brief index typedef
    using IndexType = typename Tin::IndexType;
    /// @brief harmonic order cache
    std::vector<IndexType> _harmOrd;
    /// @brief number of columns cache
    std::vector<IndexType> _cols;
};

} // namespace Cpu
} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
