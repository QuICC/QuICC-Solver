/**
 * @file WorlandSphEnergyRule.cpp
 * @brief Source of the Worland spherical energy quadrature
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void WorlandSphEnergyRule::computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);
      int jsize = 2*size;
      Internal::Array jgrid(jsize);
      Internal::Array jweights(jsize);
      LegendreRule::computeQuadrature(jgrid, jweights, jsize);

      igrid = jgrid.bottomRows(size);
      iweights.array() = jweights.bottomRows(size).array()*igrid.array().pow(2);
   }

}
}
}
