/**
 * @file WorlandCylEnergyRule.cpp
 * @brief Source of the Worland cylindrical energy quadrature
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Quadrature/WorlandCylEnergyRule.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void WorlandCylEnergyRule::computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);
      int jsize = size;
      Internal::Array jgrid(jsize);
      Internal::Array jweights(jsize);
      LegendreRule::computeQuadrature(jgrid, jweights, jsize);

      igrid = (jgrid.array() + MHD_MP(1.0))/MHD_MP(2.0);
      iweights.array() = jweights.array()*igrid.array()/MHD_MP(2.0);
   }

}
}
}
