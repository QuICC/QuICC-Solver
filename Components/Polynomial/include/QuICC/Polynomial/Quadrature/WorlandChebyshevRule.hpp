/** 
 * @file WorlandChebyshevRule.hpp
 * @brief Implementation of the Worland Chebyshev quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_WORLANDCHEBYSHEVRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_WORLANDCHEBYSHEVRULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of the Worland Chebyshev quadrature rule
    */
   class WorlandChebyshevRule
   {
      public:
         /**
          * @brief Compute the quadrature for r in [0, 1]
          */
         void computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size);

         /**
          * @brief Compute the quadrature for x = 2r^2 - 1 in [-1, 1]
          */
         void computeXQuadrature(internal::Array& igrid, internal::Array& iweights, const int size);
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_WORLANDCHEBYSHEVRULE_HPP
