/**
 * @file SphericalBesselRule.hpp
 * @brief Implementation of the spherical Bessel quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_SPHERICALBESSELRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_SPHERICALBESSELRULE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of the spherical Bessel quadrature rule
    */
   class SphericalBesselRule: public LegendreRule
   {
      public:
         /**
          * @brief Compute the quadrature
          */
         void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size);
   };

} // namespace Quadrature
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_QUADRATURE_SPHERICALBESSELRULE_HPP
