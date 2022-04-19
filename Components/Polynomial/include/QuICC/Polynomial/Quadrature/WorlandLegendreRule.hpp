/** 
 * @file WorlandLegendreRule.hpp
 * @brief Implementation of the Worland Legendre quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_WORLANDLEGENDRERULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_WORLANDLEGENDRERULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of the Worland Legendre quadrature rule
    */
   class WorlandLegendreRule: public LegendreRule
   {
      public:
         /**
          * @brief Compute the quadrature
          */
         void computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size);
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_WORLANDLEGENDRERULE_HPP
