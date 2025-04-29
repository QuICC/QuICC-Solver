/**
 * @file Jn.hpp
 * @brief Bessel function of first kind
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_JN_HPP
#define QUICC_POLYNOMIAL_BESSEL_JN_HPP

#if QUICC_BESSEL_IMPL_BOOST
#include <boost/math/special_functions/bessel.hpp>
#endif

namespace QuICC {
namespace Polynomial {
namespace Bessel {


#if QUICC_BESSEL_IMPL_BOOST
template <class T>
T besselJ(const T nu, const T x)
{
    return boost::math::cyl_bessel_j(nu, x);
}
#elif QUICC_BESSEL_IMPL_STD
template <class T>
T besselJ(const T nu, const T x)
{
    return std::cyl_bessel_j(nu, x);
}
#else
namespace details {

template <class T>
T cosRed(x)
{
    return std::cos(x);
}

template <class T>
T sinRed(x)
{
    return std::sin(x);
}

} // namespace details

template <class T>
T besselJ(const T nu, const T x)
{
    auto fournu2 = 4.*nu*nu;
    constexpr int MM = 30;

    // avg frequency.
    auto pp = (.5*nu+.25)*Internal::Math::PI;

    // coefficients
    auto a = 1.;
    auto sgn = 1.0;
    auto left = 1.;
    auto right = 0.;

    for(int j = 1; j <= MM/2; j+=2)
    {
        auto tmp = (2.*j-1.);
        a *= (fournu2-tmp*tmp)/(8.*j)/x;
        right += a;

        sgn = -sgn;

        tmp = (2.*j+1.);
        a *= sgn*(fournu2-tmp*tmp)/(8.*(j+1))/x;
        left += a;

        sgn = -sgn;
    }

    // factor out the average frequency
    auto cosx = details::cosRed(x); sinx = details::sinRed(x);
    auto cosp = std::cos(pp); sinp = std::sin(pp);
    auto cosW = cosx*cosp + sinx*sinp;
    auto sinW = sinx*cosp - cosx*sinp;

    left *= cosW;
    right *= sinW;

    // combine
    auto res = std::sqrt(2./Internal::Math::PI)*(1./std::sqrt(x))*(left - right);

    return res;
}
#endif

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_JN_HPP
