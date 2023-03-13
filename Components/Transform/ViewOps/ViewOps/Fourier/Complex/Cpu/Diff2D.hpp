/**
 * @file Diff2D.hpp
 * @brief Diff2D Complex cpu backend
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
#include "Operator/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {
/// @brief Cpu backend namespace
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Derived classes implements the differentiation in modal space
/// d^{Ofi}/d_i * d^{Ofj}/d_j + d^{Osi}/d_i * d^{Osj}/d_j
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Ofi Order of differentiation, first term, first direction
/// @tparam Ofj Order of differentiation, first term, second direction
/// @tparam Osi Order of differentiation, second term, first direction
/// @tparam Osj Order of differentiation, second term, second direction
/// @tparam Direction Fft direction
/// @tparam Treatment special treatment tag, typically of mode zero
template<class Tout, class Tin, std::size_t Ofi, std::size_t Ofj,
    std::size_t Osi, std::size_t Osj, class Direction, class Treatment = void>
class Diff2DOp : public BaseOp<Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>, Tout, Tin> {
public:
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    Diff2DOp(ScaleType scale);
    /// @brief Default constructor
    Diff2DOp() = default;
    /// @brief dtor
    ~Diff2DOp() = default;
private:
    /// @brief action implementation
    /// @param out diff modes
    /// @param in modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief Give access to base class
    friend BaseOp<Diff2DOp<Tout, Tin, Ofi, Ofj, Osi, Osj, Direction, Treatment>, Tout, Tin>;
    /// @brief Scaling factor, i.e. domain size
    ScaleType mScale{1.0};
};

} // namespace Cpu
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
