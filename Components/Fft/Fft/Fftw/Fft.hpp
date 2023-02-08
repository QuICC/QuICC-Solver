/**
 * @file Fft.hpp
 * @brief Fftw backend
 */
#pragma once

// External includes
//
#include <complex>

// Project includes
//
#include "Operator/Interface.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief This namespace provides base methods and classes for Fft backends
namespace Fft {
/// @brief This namespace provides base methods and classes for Fftw backend
namespace Fftw {

using namespace QuICC::Operator;
using namespace QuICC::Memory;

/// @brief this is the generic base class for Fftw backends
/// @tparam Tout
/// @tparam Tin
template<class Tout, class Tin>
class FftOp : public BaseOp<FftOp<Tout, Tin>, Tout, Tin>
{
public:
    /// @brief action implementation, default to no implementation
    /// @param out
    /// @param in
    void applyImpl(Tout& out, const Tin& in){throw std::logic_error("Requested Fftw Operator not implemented");};
};

/// @brief Fully dense, column major 2D
using dense2D = Attributes<DimLevelType<dense_t, dense_t>>;
/// @brief Complex  tensor, input modes view type
using Cmods_t = View<std::complex<double>, dense2D>;
/// @brief Real tensor, output phys view type
using Rphys_t = View<double, dense2D>;
/// @brief Complex tensor, output phys view type
using Cphys_t = View<std::complex<double>, dense2D>;

/// @brief Complex to Real batched Fft
/// @tparam AttIn attributes describing the input View
/// @tparam AttOut attributes describing the output View
template<class AttIn, class AttOut>
class FftOp<View<double, AttOut>, View<std::complex<double>, AttIn>> :
    public BaseOp<FftOp<View<double, AttOut>, View<std::complex<double>, AttIn>>,
        View<double, AttOut>, View<std::complex<double>, AttIn>>
{
public:
    /// @brief ctor
    FftOp();
    /// @brief dtor
    ~FftOp();
    /// @brief action implementation
    /// @param out output View
    /// @param in input View
    void applyImpl(View<double, AttOut>& out, const View<std::complex<double>, AttIn>& in);
private:
    /// @brief pointer to store the fft plan
    void* _plan{nullptr};
};

/// @brief Complex to Complex batched Fft
/// @tparam AttIn attributes describing the input View
/// @tparam AttOut attributes describing the output View
template<class AttIn, class AttOut>
class FftOp<View<std::complex<double>, AttOut>, View<std::complex<double>, AttIn>> :
    public BaseOp<FftOp<View<std::complex<double>, AttOut>, View<std::complex<double>, AttIn>>,
        View<std::complex<double>, AttOut>, View<std::complex<double>, AttIn>>
{
public:
    /// @brief ctor
    FftOp();
    /// @brief dtor
    ~FftOp();
    /// @brief action implementation
    /// @param out output View
    /// @param in input View
    void applyImpl(View<std::complex<double>, AttOut>& out, const View<std::complex<double>, AttIn>& in);
private:
    /// @brief pointer to store the fft plan
    void* _plan{nullptr};
};

} // namespace Fftw
} // namespace Fft
} // namespace QuICC
