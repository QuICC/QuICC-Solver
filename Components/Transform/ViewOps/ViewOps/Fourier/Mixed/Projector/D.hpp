/**
 * @file D.hpp
 * @brief Mixed projector diff operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Operator/Unary.hpp"
#include "Operator/Binary.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
/// @brief namespace for mixed projectors (modal to physical space)
namespace Projector {

using namespace QuICC::Operator;

/// @brief This class implements a Fourier differentiation and projection
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam FftBackend  type of FFT operator
/// @tparam DiffBackend type of operator
template<class Tout, class Tin, class FftBackend, class DiffBackend>
class DOp : public UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin> {
public:
    /// @brief type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType;
    /// @brief constructor with user defined scaling factor
    /// @param scale
    DOp(ScaleType scale);
    /// @brief default constructor
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
    std::unique_ptr<UnaryOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<BinaryOp<Tin, Tin, ScaleType>> mDiff;
    /// @brief give access to base class
    friend UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin>;
};

} // namespace Projector
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
