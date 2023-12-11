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
#include "Operator/Binary.hpp"

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
/// @tparam Treatment special treatment mask, typically of mode zero or dealiasing
template<class Tout, class Tin, std::size_t Order, class Direction,  std::uint16_t Treatment = 0>
class DiffOp : public BinaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin, typename Tout::ScalarType::value_type> {
public:
    /// @brief Type of treatment mask
    static constexpr std::uint16_t TreatmentValue = Treatment;
    /// @brief Order of differentiation
    static constexpr std::size_t OrderValue = Order;
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
    /// @param fftScaling fft scaling (inverse number of grid points)
    void applyImpl(Tout& out, const Tin& in, const ScaleType fftScaling);
    /// @brief Give access to base class
    friend BinaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin, ScaleType>;
    /// @brief Scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cuda
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
