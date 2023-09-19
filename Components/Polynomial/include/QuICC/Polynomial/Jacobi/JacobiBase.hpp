/**
 * @file JacobiBase.hpp
 * @brief Implementation of the Jacobi polynomial base
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP
#define QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Precision.hpp"
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
internal::Array P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);
   // cs(0) = precision::sqrt(MHD_MP(2.0))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));
   // possibly add a switch based on a,b to use lgamma
   cs(0) = precision::sqrt(precision::pow(MHD_MP(2.0), -a-b-MHD_MP(1.0))*precision::tgamma(a+b+MHD_MP(2.0))
         / (precision::tgamma(a+MHD_MP(1.0))*precision::tgamma(b+MHD_MP(1.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = MHD_MP(0.5) * precision::sqrt(precision::pow(MHD_MP(2.0), -a-b-MHD_MP(1.0))
         * precision::tgamma(a+b+MHD_MP(2.0))*(a+b+MHD_MP(3.0))
         / (precision::tgamma(a+MHD_MP(2.0))*precision::tgamma(b+MHD_MP(2.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -(MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))*precision::sqrt((n - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/(n + a + b));
   if (n > MHD_MP(2.0))
   {
      cs(0) *= precision::sqrt((n + a + b - MHD_MP(1.0))/(MHD_MP(2.0)*n + a + b - MHD_MP(3.0)));
   }
   cs(1) = ((MHD_MP(2.0)*n + a + b)/MHD_MP(2.0))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
   cs(2) = (a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
   cs(3) = precision::sqrt((MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*(n + a)*(n + b)));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array dP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(2.0)*precision::sqrt(MHD_MP(2.0)*(a+b))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=2 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array dP2ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);

   cs(2) = precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/(MHD_MP(2.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b)));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array dPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*precision::sqrt((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(n + a + b - MHD_MP(1.0))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
   cs(3) = precision::sqrt((n + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/((n + a)*(n + b)));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for unit normalization, recursion on P
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array dPnabPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   assert(false);

   return cs;
}

/**
 * @brief Second derivative n=0 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d2P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(8.0)*precision::sqrt((a + b)*(a + b - MHD_MP(1.0)))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d2P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (precision::sqrt(MHD_MP(3.0))/MHD_MP(2.0))*precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(1.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((n+1)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)))*precision::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
   cs(3) = precision::sqrt((n + MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(2.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative n=0 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d3P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(16.0)*precision::sqrt(MHD_MP(3.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0)))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative n=1 normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d3P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(2.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative normalizer for unit normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::unity_t, TAG>::value, bool>::type = 0
>
internal::Array d3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((n+2)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))));
   cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)));
   cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))));
   cs(3) = precision::sqrt((n + MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(3.0))));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(1.0);

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = MHD_MP(0.5);

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Polynomial normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))
         / (n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n*(n + a + b));
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0);

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array dP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(0.5)*(a + b + MHD_MP(2.0));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative n=2 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array dP2ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (MHD_MP(3.0)+a+b)*(MHD_MP(4.0)+a+b);
   cs(1) = ((MHD_MP(3.0)+a)*a-(MHD_MP(3.0)+b)*b);
   cs(2) = (MHD_MP(0.25));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array dPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(1.0));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief First derivative normalizer for natural normalization, recursion on P
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array dPnabPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = (MHD_MP(2.0)*n+a+b);
   cs(1) = n*(a-b);
   cs(2) = -n*(MHD_MP(2.0)*n+a+b);
   cs(3) = MHD_MP(2.0)*(n+a)*(n+b);

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d2P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(0.25)*(a + b)*(a + b - MHD_MP(1.0));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Second derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d2P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(1.0)));

   assert(!precision::isnan(cs.sum()));

   return cs;
}

/**
 * @brief Third derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(2.0));

   assert(!precision::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Third derivative n=0 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d3P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(1);

   cs(0) = MHD_MP(0.125)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0));

   assert(!precision::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Third derivative n=1 normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d3P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(3);

   cs(0) = (a + b + MHD_MP(2.0));
   cs(1) = (a - b);
   cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(2.0)));

   assert(!precision::isnan(cs.sum()));

   return cs;
}


/**
 * @brief Second derivative normalizer for natural normalization
 */
template <class TAG,
   typename std::enable_if<std::is_same<Quadrature::natural_t, TAG>::value, bool>::type = 0
>
internal::Array d3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{
   internal::Array cs(4);

   cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
   cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
   cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(3.0));

   assert(!precision::isnan(cs.sum()));

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
internal::MHDFloat normFact(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

} // namespace JacobiBase
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_JACOBIBASE_HPP
