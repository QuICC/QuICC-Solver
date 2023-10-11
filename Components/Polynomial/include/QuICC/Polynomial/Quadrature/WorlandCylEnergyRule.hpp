/**
 * @file WorlandCylEnergyRule.hpp
 * @brief Implementation of the Worland cylindrical energy quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_WORLANDCYLENERGYRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_WORLANDCYLENERGYRULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of the Worland cylindrical energy quadrature rule
    */
   class WorlandCylEnergyRule: public LegendreRule
   {
      public:
         /**
          * @brief Compute the quadrature
          */
         void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size);
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_WORLANDCYLENERGYRULE_HPP
