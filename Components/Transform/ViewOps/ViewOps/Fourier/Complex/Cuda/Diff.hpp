/**
 * @file Diff.hpp
 * @brief Diff Complex Cuda backend
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Operator/Unary.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
/// @brief Cuda backend namespace
namespace Cuda {

using namespace QuICC::Operator;

/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment mask, typically of mode zero or dealiasing
template<class Tout, class Tin, std::size_t Order, class Direction, std::uint16_t Treatment = 0>
class DiffOp : public UnaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin> {
public:
    /// @brief Type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    DiffOp(ScaleType scale);
    /// Default constructor
    DiffOp() = default;
    /// dtor
    ~DiffOp() = default;
private:
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief give access to base class
    friend UnaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin>;
    /// @brief scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cuda
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
