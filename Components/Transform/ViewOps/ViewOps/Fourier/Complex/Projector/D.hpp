/**
 * @file D.hpp
 * @brief Complex projector diff operator
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
/// @brief namespace for complex operators
namespace Complex {
/// @brief namespace for complex projectors (modal to physical space)
namespace Projector {

using namespace QuICC::Operator;

/// @brief This class implements a Fourier differentiation and projection
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam FftBackend  type of FFT operator
/// @tparam DiffBackend type of operator
template<class Tout, class Tin, class FftBackend, class DiffBackend>
class DOp : public BaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin> {
public:
    /// @brief Type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief Constructor with user defined scaling factor
    /// @param scale
    DOp(ScaleType scale);
    /// @brief Default constructor
    DOp();
    /// @brief dtor
    ~DOp() = default;
private:
    /// @brief action implementation that does not overwrite the input
    /// @param out differentiatied physical space coefficient
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief action implementation that might modify the input
    /// @param out differentiatied physical space coefficient
    /// @param in input modes
    void applyImpl(Tout& out, Tin& in);
    /// @brief pointer to FFT operator
    std::unique_ptr<InterfaceOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<InterfaceOp<Tin, Tin>> mDiff;
    /// @brief Give access to base class
    friend BaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin>;
};

} // namespace Projector
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
