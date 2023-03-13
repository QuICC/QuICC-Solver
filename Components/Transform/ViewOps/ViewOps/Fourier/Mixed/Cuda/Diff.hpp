/**
 * @file Diff.hpp
 * @brief Diff cpu backend
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Operator/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
/// @brief Cuda backend namespace
namespace Cuda {

using namespace QuICC::Operator;

/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment tag, typically of mode zero
template<class Tout, class Tin, std::size_t Order, class Direction, class Treatment = void>
class DiffOp : public BaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin> {
public:
    /// @brief Type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    DiffOp(ScaleType scale);
    /// @brief Default constructor
    DiffOp() = default;
    /// @brief dtor
    ~DiffOp() = default;
private:
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief Give access to base class
    friend BaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin>;
    /// @brief Scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cuda
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
