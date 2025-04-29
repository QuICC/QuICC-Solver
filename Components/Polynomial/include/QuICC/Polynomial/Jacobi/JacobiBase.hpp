/**
 * @file JacobiBase.hpp
 * @brief Implementation of the Jacobi polynomial base
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP
#define QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/traits.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"

namespace QuICC {

namespace Polynomial {

namespace Jacobi {

namespace JacobiBase {

/**
 * @brief Polynomial n=0 normalizer for unit normalization
 *
 * ref: https://dlmf.nist.gov/18.3#Px1.p1
 *
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);
   // cs(0) = Internal::Math::sqrt(MHD_MP(2.0))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));
   // possibly add a switch based on a,b to use lgamma
   cs(0) = Internal::Math::sqrt(Internal::Math::pow(MHD_MP(2.0), -a-b-MHD_MP(1.0))*Internal::Math::tgamma(a+b+MHD_MP(2.0))
         / (Internal::Math::tgamma(a+MHD_MP(1.0))*Internal::Math::tgamma(b+MHD_MP(1.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = MHD_MP(0.5) * Internal::Math::sqrt(Internal::Math::pow(MHD_MP(2.0), -a-b-MHD_MP(1.0))
         * Internal::Math::tgamma(a+b+MHD_MP(2.0))*(a+b+MHD_MP(3.0))
         / (Internal::Math::tgamma(a+MHD_MP(2.0))*Internal::Math::tgamma(b+MHD_MP(2.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -(MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))*Internal::Math::sqrt((n - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/(n + a + b));
   if (n > MHD_MP(2.0))
   {
      cs(0) *= Internal::Math::sqrt((n + a + b - MHD_MP(1.0))/(MHD_MP(2.0)*n + a + b - MHD_MP(3.0)));
   }
   cs(1) = ((MHD_MP(2.0)*n + a + b)/MHD_MP(2.0))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
   cs(2) = (a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
   cs(3) = Internal::Math::sqrt((MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*(n + a)*(n + b)));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array dP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(2.0)*Internal::Math::sqrt(MHD_MP(2.0)*(a+b))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=2 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array dP2ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);

   cs(2) = Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/(MHD_MP(2.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b)));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array dPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*Internal::Math::sqrt((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(n + a + b - MHD_MP(1.0))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
   cs(3) = Internal::Math::sqrt((n + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/((n + a)*(n + b)));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for unit normalization, recursion on P
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array dPnabPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   assert(false);

   return cs;
}

/**
 * @brief Second derivative n=0 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(8.0)*Internal::Math::sqrt((a + b)*(a + b - MHD_MP(1.0)))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (Internal::Math::sqrt(MHD_MP(3.0))/MHD_MP(2.0))*Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(1.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d2Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((n+1)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)))*Internal::Math::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
   cs(3) = Internal::Math::sqrt((n + MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(2.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative n=0 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(16.0)*Internal::Math::sqrt(MHD_MP(3.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0)))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(2.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
Internal::Array d3Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((n+2)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))));
   cs(3) = Internal::Math::sqrt((n + MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(3.0))));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(1.0);

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = MHD_MP(0.5);

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))
         / (n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n*(n + a + b));
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0);

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array dP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(0.5)*(a + b + MHD_MP(2.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=2 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array dP2ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (MHD_MP(3.0)+a+b)*(MHD_MP(4.0)+a+b);
   cs(1) = ((MHD_MP(3.0)+a)*a-(MHD_MP(3.0)+b)*b);
   cs(2) = (MHD_MP(0.25));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array dPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(1.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for natural normalization, recursion on P
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array dPnabPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = (MHD_MP(2.0)*n+a+b);
   cs(1) = n*(a-b);
   cs(2) = -n*(MHD_MP(2.0)*n+a+b);
   cs(3) = MHD_MP(2.0)*(n+a)*(n+b);

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(0.25)*(a + b)*(a + b - MHD_MP(1.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(1.0)));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d2Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(2.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Third derivative n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(1);

   cs(0) = MHD_MP(0.125)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Third derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(2.0)));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Second derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
Internal::Array d3Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{
   Internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(3.0));

   assert(!Internal::Math::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Normalization factor for weights
 *
 * Fast and accurate computation of gauss–legendre
 * and Gauss–Jacobi quadrature nodes and weights
 * Hale and Townsend
 *
 */
Internal::MHDFloat normFact(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b);

} // namespace JacobiBase
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP
