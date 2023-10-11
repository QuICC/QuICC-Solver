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
#include "Types/Internal/BasicTypes.hpp"
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
Internal::MHDFloat normFact(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
{

   Internal::MHDFloat ratio = MHD_MP(1.0);
   #ifdef QUICC_MULTPRECISION
      ratio = Internal::Math::exp(Internal::Math::lgamma(n+a+1) - Internal::Math::lgamma(n+a+b+1)
         + Internal::Math::lgamma(n+b+1) - Internal::Math::lgamma(n+1));
   #else
      if(n < 64)
      {
            ratio = std::tgamma(n+a+1) * std::tgamma(n+b+1)
               / (std::tgamma(n+a+b+1) * std::tgamma(n+1));
      }
      else
      {
         // Buehring asymptotic approximation
         auto epsilon = std::numeric_limits<Internal::MHDFloat>::epsilon();
         static constexpr int ulp = 2;
         constexpr std::size_t M = 20;
         for(std::size_t m = 1; m < M; ++m)
         {
            Internal::MHDFloat delta = MHD_MP(1.0);
            for(std::size_t i = 0; i < m; ++i)
            {
               delta *= (a + i)*(b + i)/((i+1)*(-n + i));
            }
            auto tol = epsilon * Internal::Math::abs(ratio) * ulp;
            if (Internal::Math::abs(delta) <= tol)
            {
               break;
            }
            ratio += delta;
         }
      }
   #endif
   return Internal::Math::pow(MHD_MP(2.0), a+b+1) * ratio;
}

} // namespace JacobiBase
} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC
