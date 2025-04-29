/**
 * @file Diff.hpp
 * @brief Diff Complex cpu backend
 */
#pragma once

// System includes
//
#include <memory>
#include <cstdint>

// Project includes
//
#include "Operator/Unary.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
/// @brief Cpu backend namespace
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment tag, typically of mode zero
template<class Tout, class Tin, std::size_t Order, class Direction, std::uint16_t Treatment = 0>
class DiffOp : public UnaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin> {
public:
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    DiffOp(ScaleType scale);
    /// @brief Default constructor
    DiffOp() = default;
    /// @brief dtor
    ~DiffOp() = default;
private:
    /// @brief action implementation
    /// @param out diff modes
    /// @param in modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief give access to base class
    friend UnaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin>;
    /// @brief scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
