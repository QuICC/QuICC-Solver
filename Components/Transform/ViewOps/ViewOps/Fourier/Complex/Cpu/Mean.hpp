/**
 * @file Mean.hpp
 * @brief Mean Complex cpu backend
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
namespace Complex {
/// @brief Cpu backend namespace
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Derived classes implement the extraction of the mean in modal space
/// @tparam Tout mean modes type
/// @tparam Tin input modes type
/// @tparam Direction Fft direction
template<class Tout, class Tin, class Direction>
class MeanOp : public BaseOp<MeanOp<Tout, Tin, Direction>, Tout, Tin> {
public:
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    MeanOp(ScaleType scale);
    /// @brief Default constructor
    MeanOp() = default;
    /// @brief dtor
    ~MeanOp() = default;
private:
    /// @brief action implementation
    /// @param out diff modes
    /// @param in modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief give access to base class
    friend BaseOp<MeanOp<Tout, Tin, Direction>, Tout, Tin>;
    /// @brief scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
