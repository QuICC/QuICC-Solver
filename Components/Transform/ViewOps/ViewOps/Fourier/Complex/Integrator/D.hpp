/**
 * @file D.hpp
 * @brief Complex Integrator diff operator
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
/// @brief namespace for complex operators
namespace Complex {
/// @brief namespace for complex Integrators (physical to modal space)
namespace Integrator {

using namespace QuICC::Operator;

/// @brief This class implements a Fourier differentiation and projection
/// @tparam Tout output modal space type
/// @tparam Tin input physical type
/// @tparam FftBackend  type of FFT operator
/// @tparam DiffBackend type of operator
template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2 = void>
class DOp : public UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>, Tout, Tin> {
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
    /// @param out differentiatied modes
    /// @param in input physical space coefficient
    void applyImpl(Tout& out, const Tin& in);
    /// @brief action implementation that might modify the input
    /// @param out differentiatied modes
    /// @param in input physical space coefficient
    void applyImpl(Tout& out, Tin& in);
    /// @brief pointer to FFT operator
    std::unique_ptr<UnaryOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<UnaryOp<Tout, Tout>> mDiff;
    /// @brief pointer to differentiation operator, second step
    std::unique_ptr<UnaryOp<Tout, Tout>> mDiff2;
    /// @brief Give access to base class
    friend UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>, Tout, Tin>;
};

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
