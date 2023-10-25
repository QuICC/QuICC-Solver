/**
 * @file Fft.hpp
 * @brief VkFft backend
 */
#pragma once

// External includes
//
#include <complex>

// Project includes
//
#include "Operator/Unary.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief This namespace provides base methods and classes for Fft backends
namespace Fft {
/// @brief This namespace provides base methods and classes for VkFft backend
namespace VkFft {

using namespace QuICC::Operator;
using namespace QuICC::Memory;

/// @brief this is the generic base class for VkFft backends
/// @tparam Tout
/// @tparam Tin
template<class Tout, class Tin>
class FftOp : public UnaryBaseOp<FftOp<Tout, Tin>, Tout, Tin>
{
public:
    /// @brief action implementation, default to no implementation
    /// @param out
    /// @param in
    void applyImpl(Tout& out, const Tin& in){throw std::logic_error("Requested VkFft Operator not implemented");};
};

/// @brief Complex to Real batched Fft
/// @tparam AttIn attributes describing the input View
/// @tparam AttOut attributes describing the output View
template<class AttIn, class AttOut>
class FftOp<View<double, AttOut>, View<std::complex<double>, AttIn>> :
    public UnaryBaseOp<FftOp<View<double, AttOut>, View<std::complex<double>, AttIn>>,
        View<double, AttOut>, View<std::complex<double>, AttIn>>
{
public:
    /// @brief ctor
    FftOp() = default;
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
    public UnaryBaseOp<FftOp<View<std::complex<double>, AttOut>, View<std::complex<double>, AttIn>>,
        View<std::complex<double>, AttOut>, View<std::complex<double>, AttIn>>
{
public:
    /// @brief ctor
    FftOp() = default;
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

/// @brief Real to Complex batched Fft
/// @tparam AttIn attributes describing the input View
/// @tparam AttOut attributes describing the output View
template<class AttIn, class AttOut>
class FftOp<View<std::complex<double>, AttOut>, View<double, AttIn>> :
    public UnaryBaseOp<FftOp<View<std::complex<double>, AttOut>, View<double, AttIn>>,
        View<std::complex<double>, AttOut>, View<double, AttIn>>
{
public:
    /// @brief ctor
    FftOp() = default;
    /// @brief dtor
    ~FftOp();
    /// @brief action implementation
    /// @param out output View
    /// @param in input View
    void applyImpl(View<std::complex<double>, AttOut>& out, const View<double, AttIn>& in);
private:
    /// @brief pointer to store the fft plan
    void* _plan{nullptr};
};

} // namespace VkFft
} // namespace Fft
} // namespace QuICC
