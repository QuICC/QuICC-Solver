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
#include "Profiler/Interface.hpp"

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
/// @tparam DiffBackend2 type of operator
/// second optional stage needed for Df1InvLapl2D
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
    /// @brief pointer to FFT operator
    std::unique_ptr<UnaryOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<UnaryOp<Tout, Tout>> mDiff;
    /// @brief pointer to differentiation operator, second step
    std::unique_ptr<UnaryOp<Tout, Tout>> mDiff2;
    /// @brief Give access to base class
    friend UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>, Tout, Tin>;
};

template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::DOp(ScaleType scale) : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>(scale))
{
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2 = std::make_unique<DiffBackend2>(scale);
    }
}

template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::DOp() : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>())
{
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2 = std::make_unique<DiffBackend2>();
    }
}


template<class Tout, class Tin, class FftBackend, class DiffBackend, class DiffBackend2>
void DOp<Tout, Tin, FftBackend, DiffBackend, DiffBackend2>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("DOp::applyImpl");

    // FFT
    mFft->apply(out, in);

    // differentiate in place
    mDiff->apply(out, out);

    // second optional stage needed for Df1InvLapl2D
    if constexpr(!std::is_same_v<DiffBackend2, void>)
    {
        mDiff2->apply(out, out);
    }
}

} // namespace Integrator
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
