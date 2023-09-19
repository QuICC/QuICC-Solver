/**
 * @file WorlandSphEnergyRule.hpp
 * @brief Implementation of the Worland spherical energy quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_WORLANDSPHENERGYRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_WORLANDSPHENERGYRULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of the Worland spherical energy quadrature rule
    */
   class WorlandSphEnergyRule: public LegendreRule
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

#endif // QUICC_POLYNOMIAL_QUADRATURE_WORLANDSPHENERGYRULE_HPP
