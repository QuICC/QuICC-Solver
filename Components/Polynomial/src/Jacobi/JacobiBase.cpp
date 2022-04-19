/** 
 * @file JacobiBase.cpp
 * @brief Implementation of the Jacobi polynomial base
 */

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
#include "QuICC/Polynomial/Jacobi/JacobiBase.hpp"
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/Quadrature/traits.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"


namespace QuICC {

namespace Polynomial {

namespace Jacobi {

namespace JacobiBase {

/**
 * @brief Normalization factor for weights
 *
 * Refs
 *
 * [1] Fast and accurate computation of Gauss–Legendre and Gauss–Jacobi quadrature nodes and weights
 *     Nicholas Hale and Alex Townsend 2013
 *
 */
internal::MHDFloat normFact(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
{

   internal::MHDFloat ratio = MHD_MP(1.0);
   #ifdef QUICC_MULTPRECISION
      ratio = precision::exp(precision::lgamma(n+a+1) - precision::lgamma(n+a+b+1)
         + precision::lgamma(n+b+1) - precision::lgamma(n+1));
   #else
      if(n < 64)
      {
            ratio = std::tgamma(n+a+1) * std::tgamma(n+b+1)
               / (std::tgamma(n+a+b+1) * std::tgamma(n+1));
      }
      else
      {
         // Buehring asymptotic approximation
         auto epsilon = std::numeric_limits<internal::MHDFloat>::epsilon();
         static constexpr int ulp = 2;
         constexpr std::size_t M = 20;
         for(std::size_t m = 1; m < M; ++m)
         {
            internal::MHDFloat delta = MHD_MP(1.0);
            for(std::size_t i = 0; i < m; ++i)
            {
               delta *= (a + i)*(b + i)/((i+1)*(-n + i));
            }
            auto tol = epsilon * precision::abs(ratio) * ulp;
            if (precision::abs(delta) <= tol)
            {
               break;
            }
            ratio += delta;
         }
      }
   #endif
   return precision::pow(MHD_MP(2.0), a+b+1) * ratio;
}

} // namespace JacobiBase
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC
